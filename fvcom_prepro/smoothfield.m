function field = smoothfield(fieldin,Mobj,SmoothFactor,nLoops,SmoothPts)
% Smooth a vertex-based field using averages  
%
% [field] = function smoothfield(fieldin,Mobj,SmoothFactor,nLoops,SmoothPts)  
%
% DESCRIPTION:
%    Smooth a vertex based field
%
% INPUT
%    Mobj         = Matlab mesh object
%    fieldin      = vertex based field
%    SmoothFactor = smoothing factor (0, no smoothing, 1 full smoothing)
%    nLoops       = number of smoothing iterations
%    SmoothPts    = list of vertices to smooth [optional, default = all]
%
% OUTPUT:
%    field = smoothed, vertex-based field
%
% EXAMPLE USAGE
%    Mobj.h = smoothfield(Mobj.h,Mobj,0.5,4)
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%   2014-02-25 Add verbose flag and some minor cosmetic changes.
%   
%==========================================================================

global ftbverbose;
[~, subname] = fileparts(mfilename('fullpath'));
if ftbverbose
    fprintf('\nbegin : %s\n', subname)
end

%--------------------------------------------------------------------------
% Parse input
%--------------------------------------------------------------------------

if exist('fieldin', 'var') * exist('Mobj', 'var') * exist('SmoothFactor', 'var') * exist('nLoops', 'var') == 0
	error('Arguments to smoothfield are missing')
end

if exist('SmoothPts', 'var')
	nPts = length(SmoothPts);
else
	nPts       = Mobj.nVerts;
	SmoothPts  = 1:Mobj.nVerts;
end

if ~Mobj.have_mets
	error('Cannot smooth field, need mesh metrics for smoothing, use setup_metrics.')
end

%--------------------------------------------------------------------------
% Smoothing Loops
%--------------------------------------------------------------------------

for ll = 1:nLoops
	field = fieldin;
	for ii = 1:nPts
  		i = SmoothPts(ii);
        ss = 0;
        for k = 1:Mobj.ntsn(i)
    		node = Mobj.nbsn(i,k);
    		ss = ss + field(node)/real(Mobj.ntsn(i));
        end
  		fieldin(i) = (1-SmoothFactor)*field(i) + SmoothFactor*ss;
    end
end
field = fieldin; 

if ftbverbose
    fprintf('end   : %s\n', subname)
end
