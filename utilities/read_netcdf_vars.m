function [M] = read_netcdf_vars(varargin)
%
% Function to read in listed variables from a netCDF file
%
% [M] = read_netcdf_vars(varargin)
%
% DESCRIPTION:
%   Optionally select a dimension to select a range for
%   This should work for any netCDF, but was written with FCOM output files in mind
%
% INPUT
%   Pass the variable names that you want to extract
%   [optional pair] filename, the netCDF filename
%   [optional triple] dimrange, the dimension name, the dimension range

% EXAMPLE USAGE
%   Extract variables time, x, y
%   M = read_netcdf_vars('time', 'x', 'y');
%
%   Extract variables time, x, y from filename text.nc:
%   M = read_netcdf_vars('filename', 'test.nc', 'time', 'x', 'y');
%
%   Extract variables time, x, y from filename text.nc but only for time
%   indicies 0-99:
%   M = read_netcdf_vars('filename', 'test.nc', 'dimrange', 'time', [0 100], 'time', 'x', 'y');
%
%   Extract variables time, x, y from filename text.nc but only for time
%   indicies 0-99 and siglay index of 0:
%   M = read_netcdf_vars('filename', 'test.nc', 'dimrange', 'time', [0 100], ...
%                        'dimrange', 'siglay', [0 1], 'time', 'x', 'y', 'u', 'v');
% Author(s)
%   Rory O'Hara Murray, Marine Scotland Science
%
% Revision history
%   v0 July 2013
%==========================================================================

dimrange = false;

% look for some keywords with some setting after them and remember which
% index of varargin are 'taken' in freeI.
freeI = ones(size(varargin));
subsample_num = 0;
for ii=1:1:length(varargin)
    keyword  = lower(varargin{ii});
    if length(keyword)<8, continue; end
    switch(keyword(1:8))
        case 'filename'
            netcdf_filename = varargin{ii+1};
            freeI([ii ii+1]) = 0;
        case 'dimrange'
            dimrange = true;
            subsample_num = subsample_num + 1;
            subsample_dim(subsample_num) = {varargin{ii+1}}
            subsample_ran(:,subsample_num) = varargin{ii+2};
            freeI([ii ii+1 ii+2]) = 0;
    end
end

% save all the non 'taken' values of varargin as a list
% of varables to look for and extract
varnames = varargin(find(freeI));

% if there isn't a netCDF filename defined, ask the user for one
if exist('netcdf_filename')==0
    [FileName,PathName] = uigetfile('*.nc', 'Select FVCOM netCDF file to read');
    netcdf_filename = [PathName FileName];
end

% open netCDF file for reading
ncid = netcdf.open(netcdf_filename, 'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% Get global attributes
for ii=1:ngatts
    M.gattname{ii} = netcdf.inqAttName(ncid,netcdf.getConstant('NC_GLOBAL'),ii-1);
    M.gattval{ii} = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),M.gattname{ii});
end

% get dimension lengths
for jj=1:ndims
    [dimname{jj} dimsize(jj)] = netcdf.inqDim(ncid, jj-1);
end

if dimrange
    %get ID for sub_var dimension
    for ii=1:size(dimname,2)
        for jj=1:subsample_num
            if findstr(dimname{ii},subsample_dim{jj}) sDID(jj) = ii-1; end
        end
    end
end

for ii=1:size(varnames,2) % loop through all the variables to extract
    varid(ii) = netcdf.inqVarID(ncid,varnames{ii});
    
    % Get the attributes of the variables listed        
    [varname,tmp,tmp,natts] = netcdf.inqVar(ncid, varid(ii));
    for jj=0:natts-1
        attname = netcdf.inqAttName(ncid, varid(ii), jj);
        if attname(1)=='_'
            attname2 = attname(2:end);
        else 
            attname2 = attname;
        end
        aI = strfind(attname2, '-');
        attname2(aI) = '_';

        eval(['M.' varnames{ii} '_att.' attname2 ' = netcdf.getAtt(ncid, varid(ii), attname);']);
    end

    % get info about the variable in question
	[varname xtype dimids atts] = netcdf.inqVar(ncid,varid(ii));

    % Take a subset if variable is dependent on the subsample variable
    if dimrange
        
        % work out which dimensions this variable has that should be
        % subsampled, assuming it has any...
        clear I
        for jj=1:length(dimids)
            tmp = find(dimids(jj)==sDID);
            if tmp I(jj) = tmp; end
        end
    end
        
    % 'I' only exists if this variable has dimensions that should be sub
    % sampled
    if dimrange & exist('I')
        clear dim_range
        % get all the dimension lengths and make a starts one (zeros)
        dim_range(2,:) = dimsize(dimids+1); % the sizes
        dim_range(1,:) = zeros(1, size(dim_range,2));   % the start positions, i.e. zeros
        
        % redefine the dim_starts and dim array
        for jj=1:length(I)
        	if I(jj) dim_range(:,jj) = [subsample_ran(1,I(jj)); subsample_ran(end,I(jj))-subsample_ran(1,I(jj))]; end
        end
                
        M.(varnames{ii}) = double(netcdf.getVar(ncid, varid(ii), dim_range(1,:), dim_range(2,:)));
    else
        % Do not subsample if 'I' doens't exist or we are not subsampling
        M.(varnames{ii}) = double(netcdf.getVar(ncid, varid(ii)));
    end
end

% close netCDF file
netcdf.close(ncid)

return