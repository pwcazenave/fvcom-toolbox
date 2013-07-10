function write_FVCOM_bath(Mobj,filename) 

% Write bathymetry to FVCOM format bathymetry file
%
% function write_FVCOM_bath(Mobj,filename)
%
% DESCRIPTION:
%    Generate an ascii FVCOM 3.x format bathymetry from Mesh object
%
% INPUT 
%   Mobj     = Mesh object
%   filename = FVCOM bathymetry file name
%
% OUTPUT:
%    FVCOM bathymetry file: filename
%
% EXAMPLE USAGE
%    write_FVCOM_bath(Mobj,'tst_cor.dat')   
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================
subname = 'write_FVCOM_bath';
global ftbverbose
if(ftbverbose)
  fprintf('\n'); fprintf(['begin : ' subname '\n']);
end;

%------------------------------------------------------------------------------
% Parse input arguments
%------------------------------------------------------------------------------
if(exist('Mobj')*exist('filename')==0)
	error('arguments to write_FVCOM_cor are incorrect')
end;

%------------------------------------------------------------------------------
% Dump the file
%------------------------------------------------------------------------------
if(lower(Mobj.nativeCoords(1:3)) == 'car')
	x = Mobj.x;
	y = Mobj.y;
else
	x = Mobj.lon;
	y = Mobj.lat;
end;
if(Mobj.have_bath)
	if(ftbverbose); fprintf('writing FVCOM bathymetry file %s\n',filename); end;
	fid = fopen(filename,'w');
	fprintf(fid,'Node Number = %d\n',Mobj.nVerts);
	for i=1:Mobj.nVerts
	  fprintf(fid,'%f %f %f\n',x(i),y(i),Mobj.h(i));
	end;
	fclose(fid);
else
	error('can''t write bathymetry to file, mesh object has no bathymetry')
end;

if(ftbverbose)
  fprintf(['end   : ' subname '\n'])
end;


