function Mobj = get_NAE2_forcing(Mobj, inputConf)
% Get the required parameters from NAE2 data to force FVCOM (through
% Casename_wnd.nc).
% 
% Mobj = get_NAE2_forcing(Mobj, inputConf)
% 
% DESCRIPTION:
%   Extract meteorological forcing data from Met Office NAE2 data to create
%   an FVCOM forcing file.
% 
% INPUT: 
%   Mobj - MATLAB mesh object
% 
% OUTPUT:
%   Mobj - MATLAB mesh object containing meteorological forcing data for
%   FVCOM, interpolated onto an unstructured grid.
% 
% The parameters which can be obtained from the AMM/S12 model output are:
%     - u wind component (uwnd)
%     - v wind component (vwnd)
%     - Sea level pressure (slp)
% 
% In addition to these, the momentum flux is calculated from wind data.
% KJT note: took this out from Pierre's version. Implement this?
%
% Algorithm is based on Phil Hall's 'produce_netcdf_input_data.py' script,
% adapted to handle NAE2 data by Karen Thurston
% ('produce_netcdf_input_dataR.py'. This Python script is, in turn, based
% on Jenny Brown's 'MET_INT.f' Fortran script.
% 
% Author(s)
%   Karen Thurston (National Oceanography Centre Liverpool)
% 
% Revision history:
%   2012-12-05 First version
% 
%==========================================================================

subname = 'get_NAE2_forcing';

global ftbverbose;
if(ftbverbose);
  fprintf('\n')
  fprintf(['begin : ' subname '\n'])
end

%% Where are the files we'll need?
if isunix       % Unix?
    metfname = ['/bank/jane/met/',datestr(inputConf.startDate,'YYYY'),...
        '/',lower(datestr(inputConf.startDate,'mmmYY')),'nae10R.dat'];
    comprfname = '/login/jane/NAE2/metintco.cs3x.nae2.compress.2';
    setupfname = '/work/jane/cs3x/prep/setupcs3xSGIl.uda';
%     elevfname = ['/bank/jane/cs3x/sarray.uda.',...
%         datestr(inputConf.startDate,'YYYY')];
elseif ispc     % Or Windows?
    metfname = ['\\store\bank\jane\met\',datestr(inputConf.startDate,'YYYY'),...
        '\',lower(datestr(inputConf.startDate,'mmmYY')),'nae10R.dat'];
    comprfname = '\\store\kthurs\from_Jane\metintco.cs3x.nae2.compress.2';
    setupfname = '\\store\work\jane\cs3x\prep\setupcs3xSGIl.uda';
%     elevfname = ['\\store\bank\jane\cs3x\sarray.uda.',...
%         datestr(inputConf.startDate,'YYYY')];
end

% Define the parameters of the cs3x grid
cs3x_nx = 198;   % number of grid cells in the x direction
cs3x_ny = 207;   % number of grid cells in the y direction
cs3x_lonstart = -19-(5/6); % longitude of the bottom left corner
cs3x_latstart = 40+(1/9);   % latitude of the bottom left corner
cs3x_loninc = 1/6;   % grid resolution in the x direction (degrees)
cs3x_latinc = 1/9;   % grid resolution in the y direction (degrees)
cs3x_lonfin = cs3x_lonstart+(cs3x_loninc*(cs3x_nx-1));  % longitude of the top right corner
cs3x_latfin = cs3x_latstart+(cs3x_latinc*(cs3x_ny-1));  % latitude of the top right corner
cs3x_lon = cs3x_lonstart:cs3x_loninc:cs3x_lonfin;   % array of grid lon points
cs3x_lat = cs3x_latstart:cs3x_latinc:cs3x_latfin;   % array of grid lat points
% cs3x_area = cs3x_loninc*cs3x_latinc;

% Sanity check. Does our FVCOM grid fit within the cs3x domain?
if min(Mobj.lon) < cs3x_lonstart || max(Mobj.lon) > cs3x_lonfin || ...
        min(Mobj.lat) < cs3x_latstart || max(Mobj.lat) > cs3x_latfin
    error('Your FVCOM grid is bigger than the available cs3x grid. Choose another met forcing option or crop your FVCOM grid.')
end

% Open the met compression info file
fid = fopen(comprfname);
tline = fgets(fid);
count=0;
while ischar(tline)
    count=count+1;
    compr_info{count}=tline;
    tline = fgets(fid);
end
fclose(fid);

% This section reads info from the met data compression file. The first
% line is a set of constants (NCCP, NRRP, NCCW, NRRW, NMOD/NTOTI).

% Read the first line of the compression info file
temp = sscanf(char(compr_info{1}),'%u');
NCCP = temp(1);     % Number of columns in the pressure data
NRRP = temp(2);     % Number of rows in the pressure data
NCCW = temp(3);     % Number of columns in the wind data
NRRW = temp(4);     % Number of rows in the wind data
NTOTI = temp(5);    % Also referred to as 'NMOD'

% Initialise compr_array
compr_array = nan(10,length(compr_info)-1);

for i=2:length(compr_info)
    temp = sscanf(char(compr_info{i}),'%8f');
    compr_array(:,i-1) = temp;
end

nelf = max(compr_array(1:4*NTOTI))+1;   % No idea what this is for

% Index of specified corner of source gridbox (pressure)
IBLP = compr_array(1:NTOTI)-1;    % Bottom left point
IBRP = compr_array(NTOTI+1:2*NTOTI)-1;  % Bottom right point

% Index of specified corner of source gridbox (wind)
IBLW = compr_array(2*NTOTI+1:3*NTOTI)-1;    % Bottom left point
IBRW = compr_array(3*NTOTI+1:4*NTOTI)-1;    % Bottom right point

% Weight applied to value at specified corner of source gridbox (pressure)
WTRP = compr_array(4*NTOTI+1:5*NTOTI);  % Top right point
WBRP = compr_array(5*NTOTI+1:6*NTOTI);  % Bottom right point
WTLP = compr_array(6*NTOTI+1:7*NTOTI);  % Top left point
WBLP = compr_array(7*NTOTI+1:8*NTOTI);  % Bottom left point

% Weight applied to value at specified corner of source gridbox (wind)
WTRW = compr_array(8*NTOTI+1:9*NTOTI);  % Top right point
WBRW = compr_array(9*NTOTI+1:10*NTOTI);  % Bottom right point
WTLW = compr_array(10*NTOTI+1:11*NTOTI);  % Top left point 
WBLW = compr_array(11*NTOTI+1:12*NTOTI);  % Bottom left point 

% Coefficients to rotate winds onto equatorial lat/lon grid from ELF model
COEFF1 = compr_array(12*NTOTI+1:13*NTOTI);
COEFF2 = compr_array(13*NTOTI+1:14*NTOTI); 

clear compr_info compr_array

% This section reads cs3x model control data from the setup file.
fid=fopen(setupfname,'r','b');
temp=fread(fid,inf,'*int32','b');
fclose(fid);

NRR = temp(2);  % Number of rows in full rectangle
NCC = temp(3);  % Number of columns in full rectangle
ITOT = temp(4); % Total number of points in compact arrays
IINZ = temp(5); % Number of internal Z-points

NTRNS = temp(8:NRR+7);  % Transformation matrix for compact addressing
NTRNT = temp(NRR+10:2*NRR+9);   % Transformation matrix for compact addressing
NSUM = temp(2*NRR+12:3*NRR+11); % Transformation matrix for compact addressing

clear temp

% The array section numbers (as in MET_INT.F).
hw0 = 1;
hw1 = 326;
hw2 = 100;
hw3 = 521;
hw4 = 600;
hw5 = (hw1-hw0)*(hw3-hw2);

% As far as I can tell, this section generates the cross-indexing between
% the met data and the cs3x grid.
ELFI = zeros(hw5,1);
FIEJ = zeros(hw5,1);

count = 1;
j = 1;

for iy = hw0:(hw1-1)
    for ix = hw2:(hw3-1)
        i = ix+iy*hw4;
        j = j+1;
        ELFI(count) = i-1;
        FIEJ(count) = j-1;
        count = count+1;
    end
end

if max(ELFI)+1 > nelf
    nelf = max(ELFI)+1;
end

if max(FIEJ)+1 > nelf
    nelf = max(FIEJ)+1;
end

%%
% Create an array of daily timesteps, ensuring the output time series is
% at least as long as the FVCOM model run time.
timesteps = datevec(datenum(inputConf.startDate):datenum(inputConf.endDate)+1);

% Find the number of months in the timeseries
[months,ia]=unique(timesteps(:,1:2),'rows');

for i=1:size(months,1)
    % Where are the met data files?
    if isunix       % Unix?
        metfname = ['/bank/jane/met/',num2str(timesteps(ia(i),1)),...
            '/',lower(datestr(timesteps(ia(i),:),'mmmYY')),'nae10R.dat'];
    elseif ispc     % Or Windows?
        metfname = ['\\store\bank\jane\met\',num2str(timesteps(ia(i),1)),...
            '\',lower(datestr(timesteps(ia(i),:),'mmmYY')),'nae10R.dat'];
    end
    
    PASTIT = false;
    kline = 0;
    
    % Open the met file
    fid = fopen(metfname);
    
    % calculate the number of lines per data segment
    I = fgets(fid);
    m1 = sscanf(I,'%u');
    data_seg = ceil(m1(1)/m1(2));
    
    frewind(fid);
    
    if i==1
        % Initialise the met data array
        met_temp = cell(1);
    end
    
    % Get met data from the met file and write to temporary cell array
    while PASTIT == false
        kline = kline + 2;
        I = fgets(fid);
        if I==-1
            % break the loop, it's the EOF
            break
        end
        J = fgets(fid);
        if J==-1
            % break the loop, it's the EOF
            break
        end
        
        % find the date in the met file
        datem = sscanf(J,'%u');
        datem = datenum(datem(4),datem(3),datem(2),datem(1),0,0);
        
        % Compare met file date with inputConf.startDate (is it bigger?) and
        % with inputConf.endDate (is it smaller?) If yes and yes, then we want
        % this. Write data to a temporary array.
        if (datem >= datenum(inputConf.startDate)) && ...
                (datem <= datenum(inputConf.endDate))
            % Initialise temp cell array
            temp = cell(data_seg+2,1);
            % write I and J
            temp{1} = sscanf(I,'%u');
            temp{2} = sscanf(J,'%u');
            % Get the actual data
            for m = 1:data_seg
                I = fgetl(fid);
                temp{m+2} = sscanf(I,'%10f');
            end
            met_temp{end+1}=temp;
            clear temp;
        else
            PASTIT = true;
        end
    end
    
    fclose(fid);
end

met_temp = met_temp(2:end);
    
% Not sure what these numbers are for.
JBLP = IBLP-NCCP;
JBRP = IBRP-NCCP;
JBLW = IBLW-NCCW;
JBRW = IBRW-NCCW;

elf = zeros(nelf,1);

% Temporary arrays to hold pressure/u-wind/v-wind from extracted met data
P2 = zeros(NTOTI,1);
tmpa = zeros(NTOTI,1);
tmpb = zeros(NTOTI,1);

% FirstInt = true;
kline = 1;

% Arrays to hold the pressure/u-wind/v-wind data when interpolated to the
% FVCOM grid
slp = zeros(Mobj.nVerts,size(met_temp,2)/3);
uwnd = zeros(Mobj.nElems,size(met_temp,2)/3);
vwnd = zeros(Mobj.nElems,size(met_temp,2)/3);
time = zeros(size(met_temp,2)/3,1);

% Take our forcing data, translate it to the cs3x grid, then interpolate
% onto the FVCOM grid.
for m=0:3:size(met_temp,2)-1
    for ISW = 1:3
        kline = kline+2;
        ipts =  met_temp{m+ISW}{1}(1);
        kpts =  met_temp{m+ISW}{1}(2);
        
        % Store the date
        if ISW == 1
            time((m/3)+1) = datenum(met_temp{m+ISW}{2}(4),...
                met_temp{m+ISW}{2}(3),met_temp{m+ISW}{2}(2),...
                met_temp{m+ISW}{2}(1),0,0);
        end
        
        % Preallocate 'field' for speed
        field = NaN(8,data_seg);
        
        % Get the met data we're interested in
        for n=1:data_seg
            field(:,n)=met_temp{m+ISW}{n+2};
        end
        field = reshape(field,data_seg*8,1);
        
        % if the last column of data is not full, adjust for that
        if size(met_temp{m+ISW}{end},1)<8
            too_many = 8-size(met_temp{m+ISW}{end},1);
            field = field(1:end-too_many);
        end
        
        elf(ELFI)=field(FIEJ);
        
        % Interpolate the data to the sea model grid
        if ISW == 1
            P2 = (WBLP.*elf(IBLP)')+(WBRP.*elf(IBRP)')+(WTLP.*elf(JBLP)')+...
                (WTRP.*elf(JBRP)');
        elseif ISW == 2
            tmpa = (WBLW.*elf(IBLW)')+(WBRW.*elf(IBRW)')+(WTLW.*elf(JBLW)')+...
                (WTRW.*elf(JBRW)');
        elseif ISW == 3
            tmpb = (WBLW.*elf(IBLW)')+(WBRW.*elf(IBRW)')+(WTLW.*elf(JBLW)')+...
                (WTRW.*elf(JBRW)');
        end
        
        kline = kline+data_seg;
    end
    
    % Rotate the winds to equatorial lat/lon
    tmpc = (COEFF1.*tmpa)+(COEFF2.*tmpb);
    tmpd = (COEFF1.*tmpb)-(COEFF2.*tmpa);
    
    k = 0;
    
    % Convert the data to a gridded array instead of a vector
    for j = 0:NRR-1
        I1 = NTRNT(j+1);
        I2 = NTRNS(j+1);
        for i = I1+1:I2
            k = k+1;
            PRESS(i,NRR-j) = P2(k);
            U10E(i,NRR-j) = tmpc(k);
            U10N(i,NRR-j) = tmpd(k);
        end
    end
    
    % Interpolate the pressure data onto the FVCOM grid
    slp(:,(m/3)+1) = interp2(cs3x_lon,cs3x_lat,PRESS',Mobj.lon,Mobj.lat);
    
    % Calculate the coordinates of the centre of each FVCOM grid element
    xelement = mean(Mobj.lon(Mobj.tri(:,:)),2);
    yelement = mean(Mobj.lat(Mobj.tri(:,:)),2);
    
    % Interpolate the wind data onto the FVCOM grid
    uwnd(:,(m/3)+1) = interp2(cs3x_lon,cs3x_lat,U10E',xelement,yelement);
    vwnd(:,(m/3)+1) = interp2(cs3x_lon,cs3x_lat,U10N',xelement,yelement);
end

% Convert data time to Modified Julian Day time for FVCOM
time = datevec(time);
MJD_time = greg2mjulian(time(:,1),time(:,2),time(:,3),time(:,4),time(:,5),...
    time(:,6));

% Add slp, uwnd, vwnd and time to Mobj
Mobj.Met.slp.node = slp;
Mobj.Met.uwnd.data = uwnd;
Mobj.Met.vwnd.data = vwnd;
Mobj.Met.time = MJD_time;

if(ftbverbose);
    fprintf(['end   : ' subname '\n']);
end
