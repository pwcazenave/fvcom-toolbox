function [Mobj] = add_coriolis(Mobj,cortype,fval) 

% Generate latitude used for FVCOM Coriolis file 
%
% [Mobj] = function add_coriolis(Mobj,varargin)
%
% DESCRIPTION:
%    add Coriolis parameter to Matlab mesh object
%
% INPUT 
%   Mobj                   = matlab mesh object 
%   [optional] cortype      = coriolis type 
%                             'uselatitude' (default): use Mobj.lat
%                             'constant'        
%   [optional] fval         = constant latitude for constant Coriolis
%
% OUTPUT:
%    Mobj = matlab structure containing Mesh data + fvcom Coriolis
%
% EXAMPLE USAGE
%    Mobj = add_coriolis(Mobj,'constant',41.0) 
%    Mobj = add_coriolis(Mobj) 
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================
global ftbverbose
subname = 'add_coriolis';
if(ftbverbose)
  fprintf('\n')
  fprintf(['begin : ' subname '\n'])
end;

%------------------------------------------------------------------------------
% Parse arguments 
%------------------------------------------------------------------------------
CorType = 'uselatitude';
if(exist('cortype'))
	if(cortype(1:3)=='use')
		CorType = 'uselatitude';
		if(~Mobj.have_cor)
			error('To set Coriolis with latitude, need (lon,lat) field in Mesh structure')
		end;
	else
		CorType = 'constant';
		if(~exist('fval'))
			error('must provide a constant latitude for constant coriolis option');
		end;
	end;
end;

%------------------------------------------------------------------------------
% Set Coriolis
%------------------------------------------------------------------------------
if(CorType(1:3) == 'use')
	if(ftbverbose); fprintf('setting Coriolis to latitude\n');end;
	Mobj.f = Mobj.lat;
end;

if(CorType(1:3) == 'con')
	if(ftbverbose); fprintf('setting Coriolis to constant %f\n',fval); end;
	Mobj.f = fval*ones(Mobj.nVerts,1);
end;

Mobj.have_cor = true;

if(ftbverbose)
  fprintf(['end   : ' subname '\n'])
end;




