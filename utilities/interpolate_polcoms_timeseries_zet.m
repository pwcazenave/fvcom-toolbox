function [Mobj]=interpolate_polcoms_timeseries_zet(Mobj,polcoms,relaxLine)
%% Testing for parallel matlab toolbox
wasOpened = false;
if license('test', 'Distrib_Computing_Toolbox')
    % We have the Parallel Computing Toolbox, so launch a bunch of workers.
    if matlabpool('size') == 0
        % Force pool to be local in case we have remote pools available.
        matlabpool open local
        wasOpened = true;
    end
end
%% Get info on nodes and elements involved in the nesting layers
oNodes = Mobj.relaxBC_nodes{relaxLine};
ntimes=polcoms.ntimes;
ndepths=polcoms.ndepths;
% interpolate bathymetry from polcoms onto the FVCOM nodes and elements of
% the nest layer
fdb = TriScatteredInterp(polcoms.xb(:), polcoms.yb(:), polcoms.bathy(:), 'natural');
polcoms.hb=fdb(Mobj.x(oNodes),Mobj.y(oNodes));

%%
%
% Extract distance to coast at BC points
% distance has been calculated on the nodes
distbc=Mobj.dist(oNodes);

tic
parfor di = 1:ntimes
% for di = 1:ntimes
    % Set up the interpolation objects.
    fzet = TriScatteredInterp(polcoms.bcxb(:), polcoms.bcyb(:), polcoms.zet(di,:)', 'natural');
    % Interpolate variables onto the unstructured grid.
    tempzet = fzet(Mobj.x(oNodes),Mobj.y(oNodes));
    
    % fvcom will generally have values outside polcoms domain at the
    % coast... we need to extrapolate ... use distance to coast?
    if any(isnan(tempzet))
        % split bc into two sides (assumes Boundary is surrounded by coast
        % on both sides.
        % Find max distance from coast (i.e. middle point)
        % Interpolate as a function of distance from coast
        % this can go horribly wrong and results should be checked
        [tempzet]=interpolate_near_coast(distbc,tempzet,Mobj.doExtrap);
        
    end
    
    % Interpolate single level variables onto the unstructured grid.
    fvzet(:, di) =tempzet; % surface elevation on b points
    % do depth resolving variable timeseries here
end
% for tt=1:5:ntimes
% figure(1)
% pcolor(repmat(cumsum(distbc),1,size(Mobj.siglayz, 2))',Mobj.siglayz(oNodes, :)',fvtemp(:,:,tt)');colorbar
% pause (0.2)
% end
% 
% figure(2)
% pcolor(repmat(cumsum(distbc),1,polcoms.params.n-2)',fvzb(:,:,20)',fvT(:,:,20)');colorbar
if wasOpened
    matlabpool close
end
toc

Mobj.zet=fvzet;

return


