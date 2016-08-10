function write_FVCOM_cor(Mobj,filename) 
% Write coriolis to FVCOM format coriolis file
%
% function write_FVCOM_cor(Mobj,filename)
%
% DESCRIPTION:
%    Generate an ascii FVCOM 3.x format coriolis from Mesh object
%
% INPUT 
%   Mobj     = Mesh object
%   filename = FVCOM coriolis file name
%
% OUTPUT:
%    FVCOM coriolis file: filename
%
% EXAMPLE USAGE
%    write_FVCOM_cor(Mobj, 'tst_cor.dat')   
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%    Rory O'Hara Murray (Marine Scotland Science)
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%    2014-10-07 Removed loops to speed up writing the file
%    2016-08-08 Remove tabs and clean up syntax.
%   
%==============================================================================

[~, subname] = fileparts(mfilename('fullpath'));
global ftbverbose
if ftbverbose
    fprintf('\nbegin : %s\n', subname);
end

%------------------------------------------------------------------------------
% Parse input arguments
%------------------------------------------------------------------------------
if exist('Mobj')*exist('filename') == 0
    error('arguments to write_FVCOM_cor are incorrect')
end

%------------------------------------------------------------------------------
% Dump the file
%------------------------------------------------------------------------------
if lower(Mobj.nativeCoords(1:3)) == 'car'
    x = Mobj.x;
    y = Mobj.y;
else
    x = Mobj.lon;
    y = Mobj.lat;
end

if Mobj.have_cor
    if(ftbverbose); fprintf('writing FVCOM coriolis file %s\n',filename); end;
    fid = fopen(filename,'w');
    fprintf(fid,'Node Number = %d\n',Mobj.nVerts);
    fprintf(fid, '%f %f %f\n', [x, y, Mobj.f]');
    fclose(fid);
else
    error('can''t write coriolis to file, coriolis is not setup, see add_coriolis')
end

if ftbverbose
    fprintf('end   : %s\n', subname)
end


