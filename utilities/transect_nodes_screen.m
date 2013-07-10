function [transect]  = transect_nodes_screen(Mobj,VarType)
% function to select and output index of nodes along a defined transect
% By clicking on points on the screen
%   Choose between nodes or vertices (U,V or tracers)
% [transect]  = transect_nodes_screen(Mobj,VarType)
%
% DESCRIPTION:
%    Select using ginput the set of nodes comprising a transect. It will
%    interpolate between consecutive points to the resolution of the mesh
%
% INPUT
%    Mobj = Matlab mesh object
%    VarType = node or element (for elements, one needs xc and yc from
%    FVCOM NC file. Not implemented yet.
%
%
% OUTPUT:
%    transect = Matlab structure with indices
%
% EXAMPLE USAGE
%    [transect]  = transect_nodes_screen(Mobj,'nodes')
%
% Author(s):
%    Ricardo Torres and Geoff Cowles (University of Massachusetts Dartmouth)
%
% Note:
%    Uses ginput2 which allows zoom/pan before selecting points and displays
%    clicked points realtime
%
% Revision history
%
%==============================================================================
subname = 'transect_nodes_screen';
global ftbverbose
if(ftbverbose)
    fprintf('\n')
    fprintf(['begin : ' subname '\n'])
end;
%%-----------------------------------------------------------
% Plot the mesh
%------------------------------------------------------------------------------

if(lower(Mobj.nativeCoords(1:3)) == 'car')
    x = Mobj.x;
    y = Mobj.y;
else
    x = Mobj.lon;
    y = Mobj.lat;
end;

figure
patch('Vertices',[x,y],'Faces',Mobj.tri,...
    'Cdata',Mobj.h,'edgecolor','k','facecolor','interp');
hold on;
plot(x(1316),y(1316),'ro')
% use ginput2 (which allows zooming and plots points as they are clicked) to let
% user select the boundary points
[xselect] = ginput2(true,'k+')
[npts,jnk] = size(xselect);

if(npts == 0)
    fprintf('you didn''t select any points')
    fprintf(['end   : ' subname '\n'])
    return
end;
fprintf('you selected %d points\n',npts)

% snap to the closest vertices
for i=1:npts
    [ipt(i),dist] = find_nearest_pt(xselect(i,1),xselect(i,2),Mobj);
end;

% replot domain with snapped vertices
plot(x(ipt),y(ipt),'ro');
% find nodes closes to straight lines between selected points
transect.x=[];
transect.y=[];
transect.idx=[];
for i=1:npts-1
    dI=1;
    keepgoing=1;
    while keepgoing
        dx=diff(x([ipt(i),ipt(i+1)]));
        dy=diff(y([ipt(i),ipt(i+1)]));
        datadxy=[dx,dy];
        [~,direction]=max(abs([dx,dy]));
        data=circshift([x([ipt(i),ipt(i+1)]),y([ipt(i),ipt(i+1)])],[0,direction-1]);
        datadxy= circshift(  datadxy,[0,direction-1]);
        dataXI=data(1,1):datadxy(1)/dI:data(end,1);
        dataYI = interp1(data(:,1),data(:,2),dataXI);
        dataI=circshift([dataXI(:),dataYI(:)],[0,direction-1]);
        plot(dataI(:,1),dataI(:,2),'yx');
        idx=nan*ones(1,length(dataI));
        for ii=1:length(dataI)
            [idx(ii),~] = find_nearest_pt(dataI(ii,1),dataI(ii,2),Mobj);
        end;
        % find if there are duplicates in idx
        [~,igood]=unique(idx);
        if length(igood)~=length(idx)
            % there are duplicates so no need to interpolate to higher resolution
            keepgoing=0;
        else
            % increase interpolation to dI*2;
            dI=dI*2;
            disp('Increasing resolution of transect')
        end
        idx= idx(sort(igood));
        checkelems={};
        % check if there are more than two nodes per element in the transect
        % first get all possible elements in the transect
        for ii=1:length(idx)
            checkelems{ii,1}=find(Mobj.tri(:,1)==idx(ii));
            checkelems{ii,2}=find(Mobj.tri(:,2)==idx(ii));
            checkelems{ii,3}=find(Mobj.tri(:,3)==idx(ii));
            
        end
        
        unique_elems=unique(sort(cat(1,checkelems{:})));
        cc=1;newidx=zeros(length(unique_elems),2);
        for ii=1:length(unique_elems)
            test4nodes=zeros(1,3);
            for rr=1:3
                test= find(idx==Mobj.tri(unique_elems(ii),rr));
                if ~isempty(test);
                    test4nodes(rr)=test;
                end
            end
            switch length(find(test4nodes))
                case 2
                    % we have an element with two nodes. Its a keeper
                    newidx(cc,1:2)=test4nodes(find(test4nodes));
                    cc=cc+1;
                case 3
                    disp(['Too many nodes in element ',num2str(ii)])
                    disp('Removing the middle node as locations are in order')
%                     plot(x(idx(test4nodes)),y(idx(test4nodes)),'w-')
                    newidx(cc,1:2)=test4nodes([1,3]);
                     cc=cc+1;
            end
        end
        newidx=unique(newidx(:));
        newidx=newidx(find(newidx));
        if length(newidx)==length(idx)
            disp('All nodes are linked, continue to next section')
            keepgoing=0;
            
            transect.x=[transect.x;x(idx)];
            transect.y=[transect.y;y(idx)];
            transect.idx=[transect.idx;idx(:)];
            plot(transect.x,transect.y,'wo');
            
        else
            disp('Mising linked nodes in transect, continue with increased resolution')
            keepgoing=1;
        end
    end
end
% check that no repeated values are present in indices
            transect.x=unique(transect.x);
            transect.y=unique(transect.y);
            transect.idx=unique(transect.idx);
            return
