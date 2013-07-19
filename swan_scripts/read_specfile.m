function [specstruct] = read_specfile(specfile,plotall) 
% Read the ascii-based spectral distribution output from SWAN for single station 
%
% DESCRIPTION:
%   function [specstruct] = read_specfile(specfile,plotall) 
%
% INPUT 
%   specfile = swan 40.72 spectral energy file 
%   plotall  = [optional] plot all the spectral distributions
%
% OUTPUT:
%   specstruct = matlab spectral distribution structure
%
% EXAMPLE USAGE
%   read_specfile('skg4.3.spec',true) 
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================


%------------------------------------------------------------------------------
% read in the spectral file 
%------------------------------------------------------------------------------

% make sure file exists
fid = fopen(specfile,'r');
if(fid  < 0)
  error(['file: ' specfile ' does not exist']);
end;
fclose(fid);

% determine number of spectral profiles
cmd  = ['grep LOCATION ' specfile ' | grep -v LOCATIONS | wc -l > out'];
system(cmd);
fid = fopen('out','r');
C=textscan(fid,'%d',1);
nProfiles = C{1};
fclose(fid);

% read the header
fid = fopen(specfile,'r');
lin = fgetl(fid); %header
lin = fgetl(fid); %header
lin = fgetl(fid); %header
lin = fgetl(fid); %header
lin = fgetl(fid); %header
lin = fgetl(fid); %header
lin = fgetl(fid); %header


% read the data location
C = textscan(fid, '%f %f\n',1); 
x = C{1};
y = C{2};
fprintf('probe location: %f %f\n',x,y);

% number of frequencies
C = textscan(fid, '%s %s %s %s %s',1);  %AFREQ  absolute frequencies in Hz
C = textscan(fid, '%d %s %s %s',1); 
nFreq = C{1}

% report dimensions
fprintf('number of profiles %d\n',nProfiles);
fprintf('number of freqs    %d\n',nFreq);

% allocate space for data
freq = zeros(nFreq,1);          %frequencies
time = zeros(nProfiles,1);      %time (modified julian day)
dens = zeros(nProfiles,nFreq);  %spectral density in m^2/hz
wdir = zeros(nProfiles,nFreq);  %average wave direction (Cartesian) at that freq
sprd = zeros(nProfiles,nFreq);  %directional spreading

% read the frequencies
for i=1:nFreq
  C = textscan(fid, '%f', 1);
  freq(i) = C{1};
end;

% read in stuff
C = textscan(fid, '%s', 1);
lin = fgets(fid); 
lin = fgets(fid); 
lin = fgets(fid); 
lin = fgets(fid); 
lin = fgets(fid); 
lin = fgets(fid); 
lin = fgets(fid); 
lin = fgets(fid); 
lin = fgets(fid); 
lin = fgets(fid); 
lin = fgets(fid); 

% loop over probes, reading profiles
for i=1:nProfiles
  fprintf('reading profile %d\n',i)
  C = textscan(fid, '%s %s %s %s', 1);
  C = textscan(fid, '%s %d', 1);
  for j=1:nFreq
    C = textscan(fid, '%f %f %f', 1);
    pwr = C{1}; if(pwr < 0) pwr = NaN; end; 
    dir = C{2}; if(pwr < -900) dir = NaN; end; 
    spr = C{3}; if(pwr < 0) spr = NaN; end; 
    dens(i,j) = pwr;
    wdir(i,j) = dir;
    sprd(i,j) = spr;
  end;
end;

% plot option
if(plotall)
  for i=1:2:nProfiles
    plot(1./freq,dens(i,:)); hold on;
  end;
end;
