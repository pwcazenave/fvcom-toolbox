% Take my PhD shelf model grid converted to SMS format and from it make the
% necessary files for an FVCOM model.

matlabrc
close all
clc

basename = 'ukerc_v6_manual_bathy'; % utm30n
% basename = 'ukerc_v8'; % lat/long

base = '/data/medusa/pica/models/FVCOM/shelf/';

% base = ['/users/modellers/pica/Work/pml-git/data/',basename];

% basename = 'ukerc_v1';
% base = '/data/medusa/pica/models/FVCOM/shelf/';

% basename = 'tamar_final_LRv2';
% base = '/tmp/shelf/input/configs/tamar3';


% Temporary paths when on Riqui's machine
addpath('/users/modellers/pica/Work/FVCOM/fvcom-toolbox/utilities')
addpath('/users/modellers/pica/Work/FVCOM/fvcom-toolbox/fvcom_prepro/')
addpath '/tmp/shelf/matlab/'; % for the temporary set_spectide.m file.

% Create the bathymetry. Uses my modified read_sms_mesh version to add the
% obc_nodes to the mesh object.
% Spherical
shelfMesh = read_sms_mesh(...
    '2dm',fullfile(base,'raw_data/',[basename,'_latlong.2dm']),...
    'bath',fullfile(base,'raw_data/',[basename,'.pts']),...
    'coordinate','spherical','addCoriolis',true);
% Cartesian
% shelfMesh = read_sms_mesh(...
%     '2dm',fullfile(base,'raw_data/',[basename,'.2dm']),...
%     'bath',fullfile(base,'raw_data/',[basename,'.pts']),...
%     'coordinate','cartesian','addCoriolis',true);

% Create a Coriolis file from the bathy which varies with latitude. Given
% the size of the domain, this is probably necessary.
shelfMesh = add_coriolis(shelfMesh,'uselatitude');

% Parse the open boundary nodes and add accordingly
% shelfMesh.boundaryNames = {'North','East','SouthWest'}; % Order is important!
% shelfMesh.boundaryNames = {'East','North','SouthWest'}; % Order is important!
shelfMesh.boundaryNames = {'SouthWest','East','North'}; % Order is important!
% shelfMesh.boundaryNames = {'SouthWest'}; % Order is important!
if shelfMesh.have_strings
    for i=1:size(shelfMesh.read_obc_nodes,2)
        nodeList = double(cell2mat(shelfMesh.read_obc_nodes(i)));
        shelfMesh = add_obc_nodes_list(shelfMesh,nodeList,shelfMesh.boundaryNames{i},1);
    end
end

% Create a sponge file
if shelfMesh.have_strings
    for i=1:size(shelfMesh.read_obc_nodes,2)
        nodeList = double(cell2mat(shelfMesh.read_obc_nodes(i)));
        shelfMesh = add_sponge_nodes_list(shelfMesh,nodeList,[shelfMesh.boundaryNames{i},' sponge'],20000,0.0001);
        clear nodeList
    end
end

% Create a constant roughness z0 file
shelfMesh.z0 = ones(1,shelfMesh.nElems)*0.025; % or 0.015 or 0.03 - Davies and Furnes (1980) shelf model

% Variable z0.
% z0obj = read_sms_mesh(...
%     '2dm',fullfile(base,'raw_data/',[regexprep(basename,'_bathy',''),'_seds_latlong.2dm']),...
%     'bath',fullfile(base,'raw_data/',[regexprep(basename,'_bathy',''),'_seds.pts']),...
%     'coordinate','spherical','addCoriolis',false);
% shelfMesh.d50 = z0obj.h;
% clear z0obj % don't need that now
% % Scale grain size data by 1/12 (after Soulsby, 1997) to get z0_s:
% shelfMesh.z0 = shelfMesh.d50/12;

% Estimate model time step. Supply estimated velocity (m/s) and tidal range
% (m) after the mesh object.
shelfMesh = estimate_ts(shelfMesh,3,10);
fprintf('Estimated time step:\t%.2f\n',min(shelfMesh.ts))

%% Sort up the boundary conditions. Use TPXO to get constituents.

addpath /users/modellers/pica/Work/MATLAB/toolboxes/TMD2.03/
addpath /users/modellers/pica/Work/MATLAB/toolboxes/TMD2.03/FUNCTIONS/
Model = '/users/modellers/pica/Work/MATLAB/toolboxes/TMD2.03/DATA/Model_tpxo7.2';
cd('/users/modellers/pica/Work/MATLAB/toolboxes/TMD2.03/')

addpath '/tmp/shelf/matlab/'; % for the temporary set_spectide.m file.

% Use the TMD function tmd_extract_HC.m to get harmonic constants at the
% boundary nodes.

% How many do we actually want to use at the model boundaries?
% shelfMesh.Components = {'M2','S2'};
shelfMesh.Components = {'M2','S2','N2','K2','K1','O1','P1','Q1','Mf','Mm'};
% What are their periods (seconds)?
% shelfMesh.period_obc = [44714.16, 43200];
shelfMesh.period_obc = [44714.16,43200,45570.24,43082.28,86163.84,92949.84,86637.24,96726.24,1180260,2380716];


amp=cell(1,shelfMesh.nObs);
Gph=cell(1,shelfMesh.nObs);
Depth=cell(1,shelfMesh.nObs);
constList=cell(1,shelfMesh.nObs);
for i=1:length(shelfMesh.boundaryNames)
    extractType = 'z';

    % It is possible to specify the indices of the constituents of interest
    % when calling tmd_extract_HC, but it requires knowing the order
    % they're stored in the file. Easier for me to extract the constituents
    % of interest separately. This makes it a bit slower (having to
    % interpolate all the constituents is slower than a select few), but
    % it's a bit easier to code up. 
    if shelfMesh.have_lonlat
        [amp{i},Gph{i},Depth{i},constList{i}] = tmd_extract_HC(Model,shelfMesh.lat(shelfMesh.read_obc_nodes{i}),shelfMesh.lon(shelfMesh.read_obc_nodes{i}),extractType);
    else
        % Need to convert XY to latlon.
        try % to use the handy file exchange utm2deg function
            % Make cell array of all the zones because utm2deg is a bit
            % inflexible in that regard (size of utmZones must equal size
            % of x and y).
            utmZone = {'30 U'};
            utmZones=cellfun(@(x) repmat(x,length(shelfMesh.x(shelfMesh.read_obc_nodes{i})),1),utmZone,'uni',false);
            [tmpLat,tmpLon] = utm2deg(shelfMesh.x(shelfMesh.read_obc_nodes{i}),shelfMesh.y(shelfMesh.read_obc_nodes{i}),utmZones{1});
            % Get the tidal data
            [amp{i},Gph{i},Depth{i},constList{i}] = tmd_extract_HC(Model,tmpLat,tmpLon,extractType);
        catch
            error('Can''t convert X/Y positions to lat/long, so can''t extract data from the TPXO data. Consider adding utm2deg to your PATH.')
        end
    end

    for j=1:numel(shelfMesh.Components)
        fprintf('Extracting %s... ',shelfMesh.Components{j})
        posIdx = strmatch(lower(shelfMesh.Components{j}),constList{i});
        shelfMesh.amp_obc{i}(j,:) = amp{i}(posIdx,:);
        shelfMesh.phase_obc{i}(j,:) = Gph{i}(posIdx,:); % Greenwich phase
        fprintf('done.\n')
    end
end
clear posIdx

%% Find NaNs in the boundaries

for i=1:shelfMesh.nObs
    brokenBoundary=i;
    
    nanIdx = shelfMesh.read_obc_nodes{brokenBoundary}(isnan(shelfMesh.phase_obc{brokenBoundary}(1,:)));

    nanLon = shelfMesh.lon(nanIdx);
    nanLat = shelfMesh.lat(nanIdx);

    doFig=0;
    if max(nanLon)-min(nanLon)==0
        minPos = min(nanLat);
        maxPos = max(nanLat);
        doFig=1;
    elseif max(nanLat)-min(nanLat)==0
        minPos = min(nanLon);
        maxPos = max(nanLon);
        doFig=1;
    elseif isempty(nanIdx)
        fprintf('No NaNs in %s boundary.\n',shelfMesh.boundaryNames{i})
        clear nanLon nanLat nanIdx
    else
        error('Boundaries are not linear. Won''t plot %s boundary NaNs',shelfMesh.boundaryNames{i})
    end

    if doFig
        figure
        patch('Vertices',[shelfMesh.lon,shelfMesh.lat],'Faces',shelfMesh.tri,...
                'Cdata',shelfMesh.h,'edgecolor','k','facecolor','interp');
        hold on;
        plot(shelfMesh.lon(nanIdx),shelfMesh.lat(nanIdx),'wo','LineWidth',3,'MarkerSize',12)
        plot(shelfMesh.lon(nanIdx),shelfMesh.lat(nanIdx),'ko','LineWidth',3,'MarkerSize',8)
        axis('equal','tight')
    end
end

clear doFig i j brokenBoundary

%% Write out all the required files.

base = '/tmp/shelf/';
addpath '/tmp/shelf/matlab/'; % for the temporary set_spectide.m file.

% Grid
write_FVCOM_grid(shelfMesh,fullfile(base,'input/configs/',basename,[basename,'_grd.dat']));
% Bathymetry
write_FVCOM_bath(shelfMesh,fullfile(base,'input/configs/',basename,[basename,'_dep.dat']));
% Coriolis
write_FVCOM_cor(shelfMesh,fullfile(base,'input/configs/',basename,[basename,'_cor.dat']));
% Open boundaries
write_FVCOM_obc(shelfMesh,fullfile(base,'input/configs/',basename,[basename,'_obc.dat']))
% Sponge file
write_FVCOM_sponge(shelfMesh,fullfile(base,'input/configs/',basename,[basename,'_spg.dat']))
% Bed roughness (constant or variable (see above))
write_FVCOM_z0(shelfMesh.z0,fullfile(base,'input/configs/',basename,[basename,'_z0.nc']),'bottom roughness');
% Create spectides at the boundaries
SpectralFile = fullfile(base,'input/configs/',basename,[basename,'_spectide.nc']);
% Grab number of components dynamically.
set_spectide(shelfMesh,numel(shelfMesh.Components),SpectralFile,'Shelf model spectral tidal boundary input')
% Do the temperature and salinity
write_FVCOM_tsobc(fullfile(base,'input/configs/',basename,basename),53736:1:53736+360,20) % do 2006
% example_FVCOM_tsobc(fullfile(base,'input/configs/',basename,basename),42662:1:42662+30)