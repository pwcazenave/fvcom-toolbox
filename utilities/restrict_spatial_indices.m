function FVCOM = restrict_spatial_indices(FVCOM,mask_nodes,mask_elems);

% Eliminates the FVCOM nodes and elements in the lists mask_nodes and
% mask_elems.
%
% function FVCOM = restrict_spatial_indices(FVCOM,mask_nodes,mask_elems);
%
% DESCRIPTION:
%   Loops through all variables in FVCOM and restricts spatial dimensions
%
%
% INPUT:
%   FVCOM        = Structure variable with FVCOM output data
%   mask_nodes = list of node indices to remove
%   mask_elems = list of elements indices to remove
%
% OUTPUT:
%   FVCOM with fewer nodes and elements!
%
% EXAMPLE USAGE:
%   FVCOM.temp = temperature field (node,levels,times)
%   FVCOM.u = velocity field (elements,levels,times)
%   mask_nodes = (1:300) can be the boundary/nesting nodes 
%   mask_elems = (1:400) can be the boundary/nesting elements
%    FVCOM = restrict_spatial_indices(FVCOM,mask_nodes,mask_elems);
%
% Author(s):
%   Ricardo Torres (Plymouth Marine Laboratory)

%
% Revision history:
%   2018-09-17 First version 
%
%==========================================================================


[~, subname] = fileparts(mfilename('fullpath'));

global ftbverbose
if ftbverbose
    fprintf('\nbegin : %s\n', subname)
end

vnames = fields (FVCOM);

if isfield(FVCOM,'x');
nodes = length(FVCOM.x);
elseif isfield(FVCOM,'lon')
    nodes = length(FVCOM.lon);
else
    warning('No easily identifiable variable with node dimensions positions... e.g. x/y or lon/lat and I will continue')
end

if isfield(FVCOM,'xc');
elems = length(FVCOM.xc);
elseif isfield(FVCOM,'lonc')
    elems = length(FVCOM.lonc);
else
    warning('No easily identifiable variable with element dimensions positions... e.g. xc/yc or lonc/latc and I cannot continue')
end

if exist('nodes','var')
else 
    nodes=0;
end
if exist('elems','var')
else 
    elems=0;
end

for vv=vnames'
    switch size(FVCOM.(vv{1}),1) % In FVCOM variable structure, the first dimension is always the spatial dimension if it is present
        case nodes
            disp(['Clipping variable  FVCOM.',vv{1}])
            switch ndims(FVCOM.(vv{1}))
                case 1
                    FVCOM.(vv{1})(mask_nodes)=[];
                case 2
                    FVCOM.(vv{1})(mask_nodes,:)=[];
                case 3
                    FVCOM.(vv{1})(mask_nodes,:,:)=[];
            end
        case elems
             disp(['Clipping variable  FVCOM.',vv{1}])
             switch ndims(FVCOM.(vv{1}))
                case 1
                    FVCOM.(vv{1})(mask_elems)=[];
                case 2
                    FVCOM.(vv{1})(mask_elems,:)=[];
                case 3
                    FVCOM.(vv{1})(mask_elems,:,:)=[];
            end
    end
end

            
        
