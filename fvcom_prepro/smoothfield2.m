function [field]  = smoothfield(fieldin,Mobj,nLoops,SmoothPts)

% Smooth a vertex-based field using minimum value of surrounding nodes 
%
% [field] = function smoothfield(fieldin,Mobj,nLoops,SmoothPts)  
%
% DESCRIPTION:
%    Smooth a vertex based field 
%
% INPUT
%    Mobj         = Matlab mesh object
%    fielin       = vertex based field
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
%
% Revision history
%   
%==============================================================================
subname = 'smoothfield';
%fprintf('\n')
%fprintf(['begin : ' subname '\n'])

%------------------------------------------------------------------------------
% Parse input
%------------------------------------------------------------------------------

if(exist('fieldin')*exist('Mobj')*exist('nLoops') == 0)
	error('arguments to smoothfield are missing')
end;

if(exist('SmoothPts'))
	nPts = length(SmoothPts);
else
	nPts       = Mobj.nVerts;
	SmoothPts  = 1:Mobj.nVerts;
end;

if(~Mobj.have_mets)
	error('cannot smooth field, need mesh metrics for smoothing, use setup_metrics')
end;

%------------------------------------------------------------------------------
% Smoothing Loops
%------------------------------------------------------------------------------

% initialize iteration
field = fieldin;

%iterate
for ll=1:nLoops;
	field = fieldin;
	for ii=1:nPts;
  		i = SmoothPts(ii);
  		for k=1:Mobj.ntsn(i); 
    		  node = Mobj.nbsn(i,k);
                  if(abs(field(node)) < abs(field(i))); fieldin(i) = field(node); end;
  		end;
	end;
end;
field = fieldin; 


%fprintf(['end   : ' subname '\n'])

