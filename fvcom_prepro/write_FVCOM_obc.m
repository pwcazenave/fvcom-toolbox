function write_FVCOM_obc(Mobj,filename) 

% Write FVCOM format open boundary node list file  
%
% function write_FVCOM_obc(Mobj,filename)
%
% DESCRIPTION:
%    Generate an ascii FVCOM 3.x format obc node list
%
% INPUT 
%   Mobj     = Mesh object
%   filename = FVCOM obc file name
%
% OUTPUT:
%    FVCOM obc file: filename
%
% EXAMPLE USAGE
%    write_FVCOM_obc(Mobj,'tst_spg.dat')   
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================
subname = 'write_FVCOM_obc';
global ftbverbose
if(ftbverbose)
  fprintf('\n'); fprintf(['begin : ' subname '\n']);
end;

%------------------------------------------------------------------------------
% Parse input arguments
%------------------------------------------------------------------------------
if(exist('Mobj')*exist('filename')==0)
	error('arguments to write_FVCOM_obc are incorrect')
end;

%------------------------------------------------------------------------------
% Dump the file
%------------------------------------------------------------------------------
if(ftbverbose); fprintf('writing FVCOM obc %s\n',filename); end;
fid = fopen(filename,'w');

if(Mobj.nObs==0)
	fprintf(fid,'OBC Node Number = %d\n',0);
else
	counter = 0;
	Total_Obc = sum(Mobj.nObcNodes(1:Mobj.nObs));
	fprintf(fid,'OBC Node Number = %d\n',Total_Obc);
	for i=1:Mobj.nObs
		for j=1:Mobj.nObcNodes(i)
			counter = counter + 1;
			fprintf(fid,'%d %d %d\n',counter,Mobj.obc_nodes(i,j),Mobj.obc_type(i));
		end;
	end;
end;
fclose(fid);

		
if(ftbverbose)
  fprintf(['end   : ' subname '\n'])
end;

