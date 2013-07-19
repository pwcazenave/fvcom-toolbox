function write_river_info(Mobj,filename) 

% Write information on the rivers
%
% function write_river_info(Mobj,filename)
%
% DESCRIPTION:
%    Generate an ascii format containing river information this can be used
%    to setup a RIVER name list which connects model nodes to Rivers and to the
%    NetCDF file containing those rivers
%
% INPUT 
%   Mobj     = Mesh object
%   filename = river information file
%
% OUTPUT:
%    river information file: filename
%
% EXAMPLE USAGE
%    write_FVCOM_grid(Mobj,'riv_info')   
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================
subname = 'write_river_info';
fprintf('\n'); fprintf(['begin : ' subname '\n']);

%------------------------------------------------------------------------------
% Parse input arguments
%------------------------------------------------------------------------------
if(exist('Mobj')*exist('filename')==0)
	error('arguments to write_river_info are incorrect')
end;

%------------------------------------------------------------------------------
% Dump the file
%------------------------------------------------------------------------------
fid = fopen(filename,'w');
fprintf(fid,'nRivers = %d\n',Mobj.nRivers);
if(Mobj.nRivers>0)
	for i=1:Mobj.nRivers
		fprintf(fid,'River# %d name %s #Nodes %d\n',i,char(Mobj.riv_name(i)),Mobj.nRivNodes(i));
		for j=1:Mobj.nRivNodes(i)
			fprintf(fid,'%d  ',Mobj.riv_nodes(i,j));
		end;
		fprintf('\n')
	end;
end;
		


fprintf(['end   : ' subname '\n'])

