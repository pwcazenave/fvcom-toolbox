function M = read_netcdf_vars_xy(varargin)
%
% A function to extract netcdf variables for a number of (x, y) locations
% read_netcdf_vars.m is called repeatidly
%
% help read_netcdf_vars for detailed usage info
%
% INPUT
%   Pass the variable names that you want to extract
%   [optional pair] filename, the netCDF filename
%   [optional triple] dimrange, the dimension name, the dimension range

% Author(s)
%   Rory O'Hara Murray, Marine Scotland Science
%
% Revision history
%   v0 May 2014
%==========================================================================

% defaults
nodeflag = false;
cellflag = false;

% look for some keywords with some setting after them and remember which
% index of varargin are 'taken' in freeI.
freeI = true(1,size(varargin,2));
subsample_num = 0;
for ii=1:1:length(varargin)
    keyword  = lower(varargin{ii});
    vartype = class(keyword);
    if vartype(1)~='c' | length(keyword)<10, continue; end
    switch(keyword(1:10))
        case 'nodevalues'
            nodeflag = true;
            nodevalues = varargin{ii+1};
            freeI([ii ii+1]) = false;
        case 'cellvalues'
            cellflag = true;
            cellvalues = varargin{ii+1};
            freeI([ii ii+1]) = false;
    end
end

% error checking
if cellflag & nodeflag
    if length(cellvalues) ~= length(nodevalues)
        error('The number of elements and nodes are not the same. This functions only works if they are the same :-(');
        return
    end
end

% Call read_netcdf_vars.m with all the same arguments, appart from those
% used here to define xy_nodes and cells

if nodeflag && cellflag
    for count=1:length(nodevalues);
        
        nodes=nodevalues(count);
        cells=cellvalues(count);
        M1 = read_netcdf_vars(varargin{freeI}, 'dimrange', 'node', nodes+[0 1], 'dimrange', 'nele', cells+[0 1]);
        
        if count==1
            M = M1;
        elseif count>0
            % save all the variables
            varnames = fieldnames(M1);
            for ii=4:2:length(varnames)
                dimids = M1.(varnames{ii-1}).dimids;
                if size(dimids,1)>0 && dimids(1)<=1;
                    M.(varnames{ii})(count,:,:) = M1.(varnames{ii});
                end
            end
        end
        
    end
elseif nodeflag
    for count=1:length(nodevalues);
        
        nodes=nodevalues(count);
        M1 = read_netcdf_vars(varargin{freeI}, 'dimrange', 'node', nodes+[0 1]);
        
        if count==1
            M = M1;
        elseif count>0
            % save all the variables
            varnames = fieldnames(M1);
            for ii=4:2:length(varnames)
                dimids = M1.(varnames{ii-1}).dimids;
                if dimids(1)==1;
                    M.(varnames{ii})(count,:,:) = M1.(varnames{ii});
                end
            end
        end
        
    end
elseif cellflag
    for count=1:length(cellvalues);
        
        cells=cellvalues(count);
        M1 = read_netcdf_vars(varargin{freeI}, 'dimrange', 'nele', cells+[0 1]);
        
        if count==1
            M = M1;
        elseif count>0
            % save all the variables
            varnames = fieldnames(M1);
            for ii=4:2:length(varnames)
                dimids = M1.(varnames{ii-1}).dimids;
                if dimids(1)==0;
                    M.(varnames{ii})(count,:,:) = M1.(varnames{ii});
                end
            end
        end
        
    end
end