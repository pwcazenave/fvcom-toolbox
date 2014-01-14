function [Mobj] = get_cs3x(Mobj,inputConf,SetupFile,DataFile)

% Extract boundary forcing information from NOC Operational Surge Model
% output.
%
% function get_cs3x(Mobj,inputConf,SetupFile,DataFile)
%
% DESCRIPTION:
%    Extract boundary forcing information from NOC Operational Surge
%    Model output and interpolate to FVCOM open boundary nodes.
%
% INPUT
%    Mobj          = Matlab mesh object
%    inputConf     = Structure containing start and end date and time for
%                    FVCOM run
%    SetupFile     = Location of surge model setup file
%    DataFile      = Location of surge model data file
% 
% OUTPUT:
%    Mobj.surfaceElevation = Addition to Matlab mesh object. Timeseries of 
%                            surface elevation at each open boundary point
%
% EXAMPLE USAGE
%    function get_cs3x(Mobj,inputConf,SetupFile,DataFile)
%
% Author(s):  
%    Karen Amoudry (National Oceanography Centre Liverpool)
%
% Revision history
%    2014-01-09 First version.
%   
%==============================================================================
subname = 'get_cs3x';
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

% Open the compression info file
fid=fopen(SetupFile,'r','b');
temp=fread(fid,inf,'*int32','b');
fclose(fid);

NRR = temp(2);  % Number of rows in full rectangle
NCC = temp(3);  % Number of columns in full rectangle
ITOT = temp(4); % Total number of points in compact arrays

NTRNS = temp(8:NRR+7);  % Transformation matrix for compact addressing
NTRNT = temp(NRR+10:2*NRR+9);   % Transformation matrix for compact addressing

clear temp

% Define the parameters of the cs3x grid
cs3x_lonstart = -19-(5/6); % longitude of the bottom left corner
cs3x_latstart = 40+(1/9);   % latitude of the bottom left corner
cs3x_loninc = 1/6;   % grid resolution in the x direction (degrees)
cs3x_latinc = 1/9;   % grid resolution in the y direction (degrees)
cs3x_lonfin = double(cs3x_lonstart+(cs3x_loninc*(NCC-1)));  % longitude of the top right corner
cs3x_latfin = double(cs3x_latstart+(cs3x_latinc*(NRR-1)));  % latitude of the top right corner
cs3x_lon = cs3x_lonstart:cs3x_loninc:cs3x_lonfin;   % array of grid lon points
cs3x_lat = cs3x_latstart:cs3x_latinc:cs3x_latfin;   % array of grid lat points
% cs3x_area = cs3x_loninc*cs3x_latinc;  % I'll need this if I do velocity

% Sanity check. Does our FVCOM grid fit within the cs3x domain?
if min(Mobj.lon) < cs3x_lonstart || max(Mobj.lon) > cs3x_lonfin || ...
        min(Mobj.lat) < cs3x_latstart || max(Mobj.lat) > cs3x_latfin
    error('Your FVCOM grid is bigger than the available cs3x grid. Choose another met forcing option or crop your FVCOM grid.')
end

E = zeros(NCC,NRR);     % Initialise the elevation array
% U = zeros(NCC,NRR);     % Initialise the u-velocity array
% V = zeros(NCC,NRR);     % Initialise the v-velocity array

% How many files do we need to open? One per year
timesteps = datevec(datenum(inputConf.startDate):1/24:datenum(inputConf.endDate));

% Find the number of years in the timeseries
[years,~] = unique(timesteps(:,1),'rows');
%%
for p = 1:size(years,1)
    % Construct the data filename
    dfile = [DataFile,num2str(years(p))];
    
    % Open the surge model data file
    fid = fopen(dfile,'r','b');
    
    PASTIT = false;
    
    while PASTIT == false
        
        % Read the data file
        dump=fread(fid,1,'*int32','b');   % Don't need this bit
        
        if feof(fid)
            PASTIT = true;
        else
            datem = double(fread(fid,5,'*int32','b'));  % start date of data file
            
            % Start date of data file in MJD format
            mjd4 = greg2mjulian(datem(4),datem(3),datem(2),datem(1),0,0)+(datem(5)/24);
            
            if p<size(years,1) && (mjd4 >= greg2mjulian(years(p+1),1,1,0,0,0))
                PASTIT = true;
                
                % If the date is in the range we want
            elseif (mjd4 >= inputConf.startDateMJD) && (mjd4 <= inputConf.endDateMJD)
                
                dump=fread(fid,2,'*int32','b');   % Don't need this bit
                elev = fread(fid,ITOT,'*single','b');    % read in elevation
                dump=fread(fid,2,'*int32','b');   % Don't need this bit
                u = fread(fid,ITOT,'*single','b');    % read in u velocity
                dump=fread(fid,2,'*int32','b');   % Don't need this bit
                v = fread(fid,ITOT,'*single','b');    % read in v velocity
                dump=fread(fid,1,'*int32','b');   % Don't need this bit
                
                k = 0;
                
                % Reshape data into x,y array
                for j = 0:NRR-1
                    I1 = NTRNT(j+1);
                    I2 = NTRNS(j+1);
                    for i = I1+1:I2
                        k = k+1;
                        E(i,NRR-j) = elev(k);
                    end
                end
                
                % Interpolate the elevation onto the open boundary nodes
                if ~exist('SurfaceElevation','var')
                    SurfaceElevation(:,1) = interp2(cs3x_lon,cs3x_lat,E',...
                        Mobj.lon(ObcNodes),Mobj.lat(ObcNodes));
                else
                    SurfaceElevation(:,end+1) = interp2(cs3x_lon,cs3x_lat,E',...
                        Mobj.lon(ObcNodes),Mobj.lat(ObcNodes));
                end
            elseif (mjd4 < inputConf.startDateMJD)
                dump = fread(fid,7+(3*ITOT),'*int32','b');
            elseif (mjd4 > inputConf.endDateMJD)
                PASTIT = true;
            end
        end
    end
    
    fclose(fid);
end

%% Output elevation to Mobj
Mobj.surfaceElevation = SurfaceElevation;

if(ftbverbose);
    fprintf(['end   : ' subname '\n']);
end
