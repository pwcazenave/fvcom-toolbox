% Create the necessary files for FVCOM from a given SMS input unstructured
% grid.
%
% This requires the Tide Model Driver toolbox and data from
% http://polaris.esr.org/ptm_index.html, my fork of the fvcom-toolbox
% (https://github.com/pwcazenave/fvcom-toolbox) and the OPeNDAP tools for
% MATLAB (http://www.opendap.org/ml-toolbox).
%
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
%   Karen Thurston (National Oceanography Centre, Liverpool)
%
% PWC Revision history:
%
%   2012-10-08 Made a more readily understandable version which should only
%   require changes to be made at the beginning of the file.
%   2012-11-09 Add support for interpolating NCEP data onto the FVCOM grid.
%
% KJT Revision history:
%
%   2012-12-03: Added support for reading boundary tidal elevation from
%               POLCOMS/S12/AMM operational model.
%   2013-01-14: Added support for extracting meteorological forcing data
%               from NAE2 data.
%   2013-01-16: Added support for outputting sea level air pressure forcing
%               data to NetCDF file.
%   2013-01-31: Removed a number of PWC's options that we are unlikely to
%               use at NOC or aren't fully implemented yet, for this
%               tutorial version. Consult his version of this file for
%               those options.
%   2013--7-10: Added option to read in existing sigma.dat file.


%%%------------------------------------------------------------------------
%%%                          INPUT CONFIGURATION
%%%
%%% Define the model input parameters here e.g. number of tidal components,
%%% type of boundary forcing, estimated current velocity, sponge layer
%%% coefficient and radius, forcing source etc.
%%%
%%%------------------------------------------------------------------------

% Which system am I using?
% CHANGE THESE TO YOUR OWN DIRECTORIES
if isunix       % Unix?
    basedir = '/work/kthurs/';
elseif ispc     % Or Windows?
    basedir = 'R:/';    % Insert your mapped drive to \\store\login\your_name here
end

% Base folder location - where do you want your forcing input files to end
% up?
inputConf.base = [basedir,'FVCOM/v3.1.6/test/'];
inputConf.outbase = [inputConf.base,'input/'];

% Which version of FVCOM are we using (for the forcing file formats)?
inputConf.FVCOM_version = '3.1.6';

% Case name for the model inputs and outputs
% Change this to whatever you want
inputConf.casename = 'IRS_Jan08';

% Location of various useful things
if isunix       % Unix?
    % Location of grid
    inputConf.grid = '/projectsa/unstructured_grids/GRIDSTORE/GRID-FORMAT/field5.1.grid';
    % Location of S12/AMM operational model output (needs YYYY-MM.nc on end)
    inputConf.AMM_folder = '/projectsa/ISO_Modelling/POLCOMS/OPERATIONAL/OUTPUT/NETCDF/S12/AMM.hourly.';
elseif ispc     % Or Windows?
    % Location of grid
    inputConf.grid = '\\store\projectsa\unstructured_grids\GRIDSTORE\GRID-FORMAT\field5.1.grid';
    % Location of S12/AMM operational model output (needs YYYY-MM.nc on end)
    inputConf.AMM_folder = '\\store\projectsa\ISO_Modelling\POLCOMS\OPERATIONAL\OUTPUT\NETCDF\S12\AMM.hourly.';
end

% output coordinates (FVCOM only likes cartesian at the moment)
inputConf.coordType = 'cartesian'; % 'spherical' or 'cartesian'
% input coordinates (what's my input bathy in?)
inputConf.coordInput = 'spherical'; % 'spherical' or 'cartesian'
% Input grid UTM Zone (if applicable)
inputConf.utmZone = {'30 U'};

% Sponge layer parameters
inputConf.spongeRadius = -1; % in metres, or -1 for variable
inputConf.spongeCoeff = 0.001;

% z0 value in metres
inputConf.bedRoughness = 0.025; % or 0.015, 0.025 or 0.03 - Davies and Furnes (1980) shelf model

% Estimated velocity (m/s) and tidal range (m) for time step estimate
inputConf.estVel = 3;
inputConf.estRange = 10;

% Uniform temperature and salinity values
inputConf.temperature = 10;
inputConf.salinity = 35;

% Model time ([Y,M,D,h,m,s])
inputConf.modelYear = 2008;
inputConf.startDate = [inputConf.modelYear,01,01,00,00,00];
inputConf.endDate = [inputConf.modelYear,01,06,00,00,00];

% Increment used for temperature and salinity (days)
inputConf.dtTS = 1;

% Period of wind velocity change (days)
inputConf.nDays = 7;

% Open boundary nodal forcing type
%   'z' for predicted surface elevation
%   'phase-amp' for amplitudes and phases
%   'model-output' for NOC Operational Tide Surge Model output
inputConf.obcForcing = 'model-output'; 
if strcmpi(inputConf.obcForcing, 'phase-amp')
    % Use elevation harmonics, not currents from TPXO
    inputConf.extractType = 'z'; 
    % Need to cd to TPXO directory or it doesn't work
    % (Yes, this is inelegant but it's the easiest way for now)
    here = pwd; % store the current working directory to return later
    tpxo_dir = which('TMD');    % find the TPXO directory
    tpxo_dir = tpxo_dir(1:end-5);   % remove TPXO.m
    cd(tpxo_dir)    % go to TPXO directory
    % Location of the TMD model description file
    inputConf.Model = [tpxo_dir,'DATA/Model_tpxo7.2'];
elseif strcmpi(inputConf.obcForcing, 'z')
    % Need to cd to TPXO directory or it doesn't work
    % (Yes, this is inelegant but it's the easiest way for now)
    here = pwd; % store the current working directory to return later
    tpxo_dir = which('TMD');    % find the TPXO directory
    tpxo_dir = tpxo_dir(1:end-5);   % remove TPXO.m
    cd(tpxo_dir)    % go to TPXO directory
    % Location of the TMD model description file
    inputConf.Model = [tpxo_dir,'DATA/Model_tpxo7.2'];
elseif strcmpi(inputConf.obcForcing, 'model-output')
    % Use NOC Operational Tide Surge Model output
    inputConf.extractType = 'm';
end

% How many tidal constituents do we actually want to use at the model
% boundaries? Case sensitive (M2 != m2).
% Only relevant if using TPXO.
inputConf.tidalComponents = {'M2','S2','N2','K2','K1','O1','P1','Q1'};

% Give some names to the boundaries. This must match the number of node
% strings defined in SMS. Ideally, the order of the names should match the
% order in which the boundaries were made in SMS.
inputConf.boundaryNames = {'1', '2', '3'};

% Are we doing surface forcing (NCEP + OPeNDAP, or Met Office NAE),
% or just ERA Interim winds?
inputConf.doForcing = 'NAE'; % 'NCEP', 'NAE', 'ERA', 'uniform' or any
% other string to skip forcing files

%%%------------------------------------------------------------------------
%%%                      END OF INPUT CONFIGURATION
%%%------------------------------------------------------------------------

%% Prepare the data

% Convert times to Modified Julian Date
inputConf.startDateMJD = greg2mjulian(inputConf.startDate(1),inputConf.startDate(2),inputConf.startDate(3),inputConf.startDate(4),inputConf.startDate(5),inputConf.startDate(6));
inputConf.endDateMJD = greg2mjulian(inputConf.endDate(1),inputConf.endDate(2),inputConf.endDate(3),inputConf.endDate(4),inputConf.endDate(5),inputConf.endDate(6));
inputConf.inputTimeTS = inputConf.startDateMJD:inputConf.dtTS:inputConf.endDateMJD;

% Read the input mesh and bathymetry. Also creates the data necessary for
% the Coriolis correction in FVCOM.
Mobj = read_grid_mesh('grid',inputConf.grid,...
    'coordinate',inputConf.coordType,'in_coord',inputConf.coordInput,...
    'project',true,'zone',inputConf.utmZone,'addCoriolis',true);

% Parse the open boundary nodes and add accordingly
% Add the sponge nodes
for i=1:size(Mobj.read_obc_nodes,2)
    nodeList = double(cell2mat(Mobj.read_obc_nodes(i)));
    Mobj = add_obc_nodes_list(Mobj,nodeList,inputConf.boundaryNames{i},1);
    if inputConf.spongeRadius < 0    % if we want a variable sponge radius
        if i==1
            % Create an array to store the radii
            Mobj.sponge_rad = zeros(size(Mobj.sponge_nodes));
        end
        % calculate the sponge radius
        spongeRadius = calc_sponge_radius(Mobj,nodeList);
        % Add the sponge nodes to the list
        Mobj = add_sponge_nodes_list(Mobj,nodeList,...
            [inputConf.boundaryNames{i},' sponge'],spongeRadius,...
            inputConf.spongeCoeff);
    else
        Mobj = add_sponge_nodes_list(Mobj,nodeList,...
            [inputConf.boundaryNames{i},' sponge'],inputConf.spongeRadius,...
            inputConf.spongeCoeff);
    end
    clear nodeList
end

clear i

% Get the sigma depths in order to interpolate from the POLCOMS depths
if exist(fullfile(inputConf.outbase, 'sigma.dat'),'file')
    % If the sigma.dat file exists, read it
    Mobj = read_sigma(Mobj, fullfile(inputConf.base, 'input/sigma.dat'));
else
    % If we can't find the sigma.dat file, print an error message and
    % finish
    error(['sigma.dat not found. Please put your sigma.dat file into ',...
        fullfile(inputConf.outbase),' and try again.'])
end

% Do the bed roughness
Mobj.z0 = ones(1,Mobj.nElems)*inputConf.bedRoughness;

% Estimate model time step. Supply estimated velocity (m/s) and tidal range
% (m) after the mesh object.
Mobj = estimate_ts(Mobj,inputConf.estVel,inputConf.estRange);
fprintf('Estimated time step:\t%.2f\n',min(Mobj.ts))

% Extract tides from S12/AMM operational model for boundary forcing
if strcmpi(inputConf.obcForcing, 'model-output')
    % Sanity check. If model date is before 20071101:000000, can't use
    % operational model output
    if datenum(inputConf.startDate) < datenum(2007,11,1,0,0,0)
        error('Your start time pre-dates NOC Operational Tide Surge Model output availability. Check your start date, or use another boundary forcing option.')
    end
    Mobj = get_AMM(Mobj,inputConf.startDate,inputConf.endDate,inputConf.AMM_folder);
    
    % Create a time series in MATLAB datenum format with hourly inputs
    inputConf.JulianTime = datenum(inputConf.startDateMJD):1/24:datenum(inputConf.endDateMJD);
end

if strcmpi(inputConf.obcForcing,'phase-amp')    
    % Boundary conditions from TPXO (for spectral tides or predicted
	% surface elevations)

    % Put the input list into the mesh object.
    Mobj.Components = inputConf.tidalComponents;
    
    % Set up the tidal struct. This contains the relevant information for up to
    % eight constituents, ordered as period, beta love number and equilibrium
    % amplitude.
    %                    period    beta    eq. amp.
    %                      (s)      (?)       (m)
    tideComponents.M2 = [44714.16, 0.693, 0.242334];
    tideComponents.S2 = [43200.00, 0.693, 0.112841];
    tideComponents.N2 = [45570.24, 0.693, 0.046398];
    tideComponents.K2 = [43082.28, 0.693, 0.030704];
    tideComponents.K1 = [86163.84, 0.736, 0.141565];
    tideComponents.O1 = [92949.84, 0.695, 0.100514];
    tideComponents.P1 = [86637.24, 0.706, 0.046843];
    tideComponents.Q1 = [96726.24, 0.695, 0.019256];
    %tideComponents.Mf = [1180260,  ?????, 0.041742];
    %tideComponents.Mm = [2380716,  ?????, 0.022026];
    %tideComponents.Ssa = [15778980, ????, 0.019446];

    
    % Extract the values for each tidal component into Mobj.period_obc,
    % Mobj.beta_love and Mobj.equilibrium_amp.
    for c=1:size(Mobj.Components,2)
        Mobj.period_obc(c) = tideComponents.(Mobj.Components{c})(1);
        Mobj.beta_love(c) = tideComponents.(Mobj.Components{c})(2);
        Mobj.equilibrium_amp(c) = tideComponents.(Mobj.Components{c})(3);
    end
    clear c
    %
    if strcmpi(inputConf.obcForcing, 'phase-amp')
        % Provide amplitude and phase data for the boundary nodes. Use the TMD
        % function tmd_extract_HC.m to get harmonic constants at the boundary
        % nodes.
        
        amp=cell(1,Mobj.nObs);
        Gph=cell(1,Mobj.nObs);
        Depth=cell(1,Mobj.nObs);
        constList=cell(1,Mobj.nObs);
        for i=1:length(inputConf.boundaryNames)
            
            % It is possible to specify the indices of the constituents of interest
            % when calling tmd_extract_HC, but it requires knowing the order
            % they're stored in the file. Easier for me to extract the constituents
            % of interest separately. This makes it a bit slower (having to
            % interpolate all the constituents is slower than a select few), but
            % it's a bit easier to code up.
            if Mobj.have_lonlat
            [amp{i},Gph{i},Depth{i},constList{i}] = tmd_extract_HC(inputConf.Model,Mobj.lat(Mobj.read_obc_nodes{i}),Mobj.lon(Mobj.read_obc_nodes{i}),inputConf.extractType);
            else
                % Need to convert XY to latlon.
                try % to use the handy file exchange utm2deg function
                    % Make cell array of all the zones because utm2deg is a bit
                    % inflexible in that regard (size of utmZones must equal size
                    % of x and y).
                    % This is somewhat redundant now that the lat/long is added
                    % when generating the Coriolis values, but it's still
                    % worthwhile keeping it here just in case. No harm should
                    % come of it being here anyway.
                    utmZones=cellfun(@(x) repmat(x,length(Mobj.x(Mobj.read_obc_nodes{i})),1),inputConf.utmZone,'uni',false);
                    [tmpLat,tmpLon] = utm2deg(Mobj.x(Mobj.read_obc_nodes{i}),Mobj.y(Mobj.read_obc_nodes{i}),utmZones{1});
                    % Get the tidal data
                    [amp{i},Gph{i},Depth{i},constList{i}] = tmd_extract_HC(inputConf.Model,tmpLat,tmpLon,inputConf.extractType);
                catch %#ok<CTCH>
                    error('Can''t convert X/Y positions to lat/long, so can''t extract data from the TPXO data. Consider adding utm2deg to your PATH.')
                end
            end
            
            for j=1:numel(Mobj.Components)
                fprintf('Extracting %s... ',Mobj.Components{j})
                posIdx = strmatch(lower(Mobj.Components{j}),constList{i}); %#ok<MATCH2>
                Mobj.amp_obc{i}(j,:) = amp{i}(posIdx,:);
                Mobj.phase_obc{i}(j,:) = Gph{i}(posIdx,:); % Greenwich phase
                fprintf('done.\n')
            end
        end
        clear posIdx amp Gph Depth constList
        
        % Find NaNs in the boundaries
        for i=1:Mobj.nObs
            brokenBoundary=i;
            
            nanIdx = Mobj.read_obc_nodes{brokenBoundary}(isnan(Mobj.phase_obc{brokenBoundary}(1,:)));
            
            nanLon = Mobj.lon(nanIdx);
            nanLat = Mobj.lat(nanIdx);
            
            inputConf.doFig=0;
            if max(nanLon)-min(nanLon)==0
                minPos = min(nanLat);
                maxPos = max(nanLat);
                inputConf.doFig=1;
            elseif max(nanLat)-min(nanLat)==0
                minPos = min(nanLon);
                maxPos = max(nanLon);
                inputConf.doFig=1;
            elseif isempty(nanIdx)
                fprintf('No NaNs in %s boundary.\n',inputConf.boundaryNames{i})
                clear nanLon nanLat nanIdx
            else
                error('Boundaries are not linear. Won''t plot %s boundary NaNs',inputConf.boundaryNames{i})
            end
            
            if inputConf.doFig
                figure
                patch('Vertices',[Mobj.lon,Mobj.lat],'Faces',Mobj.tri,...
                    'Cdata',Mobj.h,'edgecolor','k','facecolor','interp');
                hold on;
                plot(Mobj.lon(nanIdx),Mobj.lat(nanIdx),'wo','LineWidth',3,'MarkerSize',12)
                plot(Mobj.lon(nanIdx),Mobj.lat(nanIdx),'ko','LineWidth',3,'MarkerSize',8)
                axis('equal','tight')
            end
        end
        
        clear doFig i j brokenBoundary

    elseif strcmpi(inputConf.obcForcing, 'z')
        % Use tmd_tide_pred to predict surface elevations for a given time
        % range.
        
        % Add the tidal components to the Mobj.
        Mobj.Components = inputConf.tidalComponents;
        
        % Create a time series in MATLAB datenum format with ten minute
        % inputs
        inputConf.inputTimeTideZ = datenum(inputConf.startDate):1/144:datenum(inputConf.endDate);
        % Also do Modified Julian Day for the output to NetCDF
        inputConf.inputTimeTideZMJD = datenum(inputConf.startDateMJD):1/144:datenum(inputConf.endDateMJD);
        % Get the indices to use the tidal constituents defined in
        % inputConf.tidalComponents for TPXO (which requires a numerical
        % array of the constituents to be used). The order of the TPXO
        % constituents is M2, S2, N2, K2, K1, O1, P1, Q1, MF, MM, M4, MS4,
        % MN4.
        tpxoConsts = {'M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1', ...
            'MF', 'MM', 'M4', 'MS4', 'MN4'};
        tIndUse = nan(length(Mobj.Components), 1);
        tInd = 1:length(tpxoConsts);
        for i=1:length(Mobj.Components)
            tPos = tInd(strcmp(Mobj.Components{i}, tpxoConsts));
            if ~isempty(tPos)
                tIndUse(i) = tPos;
            else
                warning('Supplied constituent (%s) is not present in the TPXO data', Mobj.Components{i}) %#ok<WNTAG>
            end
        end
        % Tidy up a bit
        clear c tpxoConsts tPos tInd

        % We can't just use tmd_tide_pred to do all the surface elevations
        % at once. Instead, the useful approaches are:
        %
        %   1. Time series at a single location
        %   2. Map of a given time step at all locations
        %
        % Since I'm likely to have many more time steps than locations,
        % it's probably best to do the time series at each location than
        % all the locations and a single time step.
        %
        % The order of the surface elevations in Mobj.surfaceElevation
        % should reflect the order of the open boundary node IDs as FVCOM
        % assumes they just map directly. So, rather than iterate through
        % each position, we need to get the position based on the list of
        % node IDs (Mobj.obc_nodes, without the zeros and in order
        % of each boundary).
        tmpObcNodes = Mobj.obc_nodes';
        ObcNodes = tmpObcNodes(tmpObcNodes~=0)';
        clear tmpObcNodes
%         obc_lat = Mobj.lat(Mobj.obc_nodes(Mobj.obc_nodes~=0));
%         obc_lon = Mobj.lon(Mobj.obc_nodes(Mobj.obc_nodes~=0));
        Mobj.surfaceElevation = nan(size(ObcNodes,2), size(inputConf.inputTimeTideZ,2));
        for i=1:size(ObcNodes,2)
            % Get the current location (from the node ID)
            currLon = Mobj.lon(ObcNodes(i));
            currLat = Mobj.lat(ObcNodes(i));
            fprintf('Position %i of %i (%.3f %.3f)... \n', i, size(ObcNodes,2), currLon, currLat)
            [Mobj.surfaceElevation(i,:), ~] = tmd_tide_pred(inputConf.Model, inputConf.inputTimeTideZ, currLat, currLon, 'z', tIndUse);
        end
        % Tidy up a bit
        clear tIndUse obc_lat obc_lon ObcNodes currLon currLat

        fprintf('done.\n')
    else
        % No idea what has been supplied
        error('Unrecognised open boundary node forcing type. Choose ''phase-amp'' or ''z'' for inputConf.obcForcing.')
    end
end

% Extract met forcing data from Met Office NAE model

if strcmpi(inputConf.doForcing, 'NAE')
    % Sanity check. If model date is before 20070501:000000, use the old
    % NAE grid
    if datenum(inputConf.startDate) < datenum(2007,5,1,0,0,0)
        % OLD NAE
        error('Met forcing data only available from 1 May 2007 onwards at the moment.\n')
    elseif datenum(inputConf.startDate) >= datenum(2007,5,1,0,0,0)
        % NAE2
        fprintf('Extracting meteorological forcing data.\n')
        Mobj = get_NAE2_forcing(Mobj,inputConf);
    end
end

%% Write out all the required files.
% Make the output directory if it doesn't exist
if exist(inputConf.outbase, 'dir')~=7
    mkdir(inputConf.outbase)
end

% Grid
write_FVCOM_grid(Mobj,fullfile(inputConf.outbase,[inputConf.casename,'_grd.dat']));

% Bathymetry
write_FVCOM_bath(Mobj,fullfile(inputConf.outbase,[inputConf.casename,'_dep.dat']));

% Coriolis
write_FVCOM_cor(Mobj,fullfile(inputConf.outbase,[inputConf.casename,'_cor.dat']));

% Open boundaries
write_FVCOM_obc(Mobj,fullfile(inputConf.outbase,[inputConf.casename,'_obc.dat']))

% Sponge file
write_FVCOM_sponge(Mobj,fullfile(inputConf.outbase,[inputConf.casename,'_spg.dat']))

% Bed roughness (constant or variable (see above))
write_FVCOM_z0(Mobj.z0,fullfile(inputConf.outbase,[inputConf.casename,'_z0.nc']),'bottom roughness');

% Tides (either spectral or surface elevation)
if strcmpi(inputConf.obcForcing, 'phase-amp')
    % Harmonic constituents from TPXO.
    SpectralFile = fullfile(inputConf.outbase,[inputConf.casename,'_spectide.nc']);
    set_spectide(Mobj,numel(Mobj.Components),SpectralFile,'TPXO spectral tidal boundary input')
elseif strcmpi(inputConf.obcForcing, 'z')
    % TPXO predicted surface elevation.
    ElevationFile = fullfile(inputConf.outbase,[inputConf.casename,'_elevtide.nc']);
    write_FVCOM_elevtide(Mobj,inputConf.inputTimeTideZMJD,ElevationFile,'Shelf model surface elevation boundary input')
elseif strcmpi(inputConf.obcForcing,'model-output')
    % Do surface elevation instead
    ElevationFile = fullfile(inputConf.outbase,[inputConf.casename,'_tides.nc']);
    write_FVCOM_elevtide(Mobj,inputConf.JulianTime,ElevationFile,'Shelf model surface elevation boundary input')
else
    error('Unrecognised open boundary node forcing type. Choose ''phase-amp'', or ''z'' for obc_type.')
end

% Do the temperature and salinity
fprintf('Writing temperature and salinity file.\n')
write_FVCOM_tsobc(fullfile(inputConf.outbase,inputConf.casename),...
    inputConf.inputTimeTS,size(Mobj.siglayz,2),inputConf.temperature,...
    inputConf.salinity)

% Write out the interpolated met forcing data to the NetCDF file
if strcmpi(inputConf.doForcing, 'NAE')
    fprintf('Writing meteorological forcing file.\n')
    
    metBase = fullfile(inputConf.outbase,inputConf.casename);
    
    write_FVCOM_forcing(Mobj,metBase,Mobj.Met,...
        [inputConf.doForcing, ' atmospheric forcing data'],...
        inputConf.FVCOM_version);
end

% Return to starting directory (if used TPXO)
if strcmpi(inputConf.obcForcing, 'phase-amp') || strcmpi(inputConf.obcForcing, 'z')
    cd(here)
end

fprintf('All done!\n')