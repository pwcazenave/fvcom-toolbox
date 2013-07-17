function [Mobj] = get_AMM(Mobj,StartDate,EndDate,ModelFolder)

% Extract boundary forcing information from NOC Operational Tide Surge
% Model output.
%
% function get_AMM(Mobj,StartDate,EndDate,ModelFolder)
%
% DESCRIPTION:
%    Extract boundary forcing information from NOC Operational Tide Surge
%    Model output and interpolate to FVCOM open boundary nodes.
%
% INPUT
%    Mobj          = Matlab mesh object
%    StartDate     = Start date and time for FVCOM run
%    EndDate       = End date and time for FVCOM run
%    ModelFolder   = Location of AMM/S12 hourly outputs
% 
% OUTPUT:
%    Mobj.surfaceElevation = Addition to Matlab mesh object
%
% EXAMPLE USAGE
%    function get_AMM(Mobj,StartDate,EndDate,ModelFolder)
%
% Author(s):  
%    Karen Thurston (National Oceanography Centre Liverpool)
%
% Revision history
%    2012-12-04 First version.
%   
%==============================================================================
subname = 'get_AMM';
global ftbverbose;
if(ftbverbose);
  fprintf('\n')
  fprintf(['begin : ' subname '\n'])
end;

%% Put the open boundary nodes into a single array for convenience
count = 0;
ObcNodes = nan(1,sum(Mobj.nObcNodes));
for ob=1:Mobj.nObs
	nObcs = Mobj.nObcNodes(ob);
    for j=1:nObcs
		count = count + 1;
		ObcNodes(count) = Mobj.obc_nodes(ob,j);  % set the open boundary nodes
    end
end

% Create an array of hourly timesteps
timesteps = datevec(datenum(StartDate):1/24:datenum(EndDate));

% Initialise an array for the sea surface elevation
SurfaceElevation = nan(count,size(timesteps,1));

% For each timestep, find the appropriate AMM/S12 file and extract the
% surface elevation
for i=1:size(timesteps,1)
    % Create AMM/S12 output filename from date
    % First, accommodate Operational Model idiosyncracy about output times
    if timesteps(i,3)==1 && sum(timesteps(i,4:6)) == 0
        tempdate = datevec(addtodate(datenum(timesteps(i,:)),-1,'month'));
        AMMfile=[ModelFolder,num2str(tempdate(1)),'-',num2str(tempdate(2),'%02i'),'.nc'];
    else
        AMMfile=[ModelFolder,num2str(timesteps(i,1)),'-',num2str(timesteps(i,2),'%02i'),'.nc'];
    end
    
    % Convert FVCOM timestep into AMM/S12 output (seconds since
    % 20071101:000000)
    AMM_time = etime(timesteps(i,:),[2007,11,01,0,0,0]);
    
    % Load the timeseries from the AMM/S12 output file
    AMM_timeseries = ncread(AMMfile,'time');
    
    % Locate the appropriate timestep number
    AMM_timestep = find(AMM_timeseries==AMM_time);

    % Load the sea surface elevation ouptut, lat and lon
    AMM_elev = ncread(AMMfile,'zet',[1 1 AMM_timestep],[Inf Inf 1])';
    AMM_lat = ncread(AMMfile,'lat');
    AMM_lon = ncread(AMMfile,'lon');
    
    % Interpolate the sea surface elevation output to the open boundary
    % nodes
    [X,Y]=meshgrid(AMM_lon,AMM_lat);
    SurfaceElevation(:,i) = interp2(X,Y,AMM_elev,Mobj.lon(ObcNodes),...
        Mobj.lat(ObcNodes));    
end

Mobj.surfaceElevation = SurfaceElevation;

if(ftbverbose);
    fprintf(['end   : ' subname '\n']);
end
