function write_FVCOM_grid(Mobj,filename) 

% Write grid and connectivity to FVCOM format grid file
%
% function write_FVCOM_grid(Mobj,filename)
%
% DESCRIPTION:
%    Generate an ascii FVCOM 3.x format gridfile from Mesh object
%
% INPUT 
%   Mobj     = Mesh object
%   filename = FVCOM grid file name
%
% OUTPUT:
%    FVCOM grid file: filename
%
% EXAMPLE USAGE
%    write_FVCOM_grid(Mobj,'tst_grd.dat')   
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%    Rory O'Hara Murray (Marine Scotland Science)
%
% Revision history
%    2014-10-07 Removed loops to speed up writing the file
%==============================================================================
subname = 'write_FVCOM_grid';
global ftbverbose
if(ftbverbose);
  fprintf('\n'); fprintf(['begin : ' subname '\n']);
end;

%------------------------------------------------------------------------------
% Parse input arguments
%------------------------------------------------------------------------------
if(exist('Mobj')*exist('filename')==0)
	error('arguments to write_FVCOM_grid are incorrect')
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
if(ftbverbose);  fprintf('writing FVCOM gridfile %s\n',filename); end;
fid = fopen(filename,'w');
fprintf(fid,'Node Number = %d\n',Mobj.nVerts);
fprintf(fid,'Cell Number = %d\n',Mobj.nElems);
% for i=1:Mobj.nElems
%   fprintf(fid,'%d %d %d %d %d\n',i,Mobj.tri(i,1:3),i);
% end;
% for i=1:Mobj.nVerts
%   fprintf(fid,'%d %f %f %f\n',i,x(i),y(i),Mobj.h(i));
% end;
fprintf(fid, '%d %d %d %d %d\n', [[1:Mobj.nElems]' Mobj.tri [1:Mobj.nElems]']');
fprintf(fid, '%d %f %f %f\n', [[1:Mobj.nVerts]' x y Mobj.h]');
fclose(fid);

if(ftbverbose);
  fprintf(['end   : ' subname '\n'])
end;

