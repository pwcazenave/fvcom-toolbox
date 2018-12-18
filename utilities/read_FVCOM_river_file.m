function [Riv]=read_FVCOM_river_file(river_root,Mobj)
% function [Riv]=read_river_file(river_root,Mobj)
% Reads FVCOM river files 
%
% function [Riv]=read_river_file(river_root,Mobj)
%
% DESCRIPTION:
%    Reads both *.nc and *nml files associated with FVCOM river input files
%
% INPUT
%   river_root = full address of river filename without extension 
%   Mobj       = needs bathymetry, nodes, elements and triangulation table. 
%
% OUTPUT:
%    Riv       = structure variable with all variables in the file 
%
% EXAMPLE USAGE
%    [Riv]=read_river_file('~/tapas_v0_riv_new',Mobj)
%
% Author(s):
%    Ricardo Torres (Plymouth Marine Laboratory) 
%
% Revision history
%
%   2018-07-22 First version.
%
%==============================================================================
[~, subname] = fileparts(mfilename('fullpath'));
global ftbverbose
if ftbverbose
    fprintf('\nbegin : %s \n', subname)
end

% Check if they are in the same order... (a couple of manual test suggest
ncRivers = [river_root '.nc'];
fidnml = fopen ([river_root '.nml']);
% get variables from netcdf file
info = ncinfo(ncRivers);
% extract all variables
varnames={};
for nn =1:length(info.Variables)
    varnames{nn}=info.Variables(nn).Name;
end
Riv=[];nn=1;
for vv=varnames
Riv.(vv{1}).data = ncread(ncRivers,vv{1});
if isfloat(Riv.(vv{1}).data)
    Riv.(vv{1}).data =double(Riv.(vv{1}).data);
    % extract units
    for rr = 1:length(info.Variables(nn).Attributes)
    if startsWith(info.Variables(nn).Attributes(rr).Name ,'units')
    Riv.(vv{1}).units =     info.Variables(nn).Attributes(rr).Value
    end 
    end
end
nn=nn+1;
end
pos =0;

while ~feof(fidnml)
    pos = pos+1;
    lines=cell(3,1);
    for dd=1:3
        lines{dd} = fgetl(fidnml);
        %     fprintf(fout,'%s\n',line)
    end
    % keep name
    [river_name]= extractAfter(lines{2},'= ');
    Riv.river_name(pos,1)={river_name(1:end-1)};
    % read grid location
    line = fgetl(fidnml);
    % extract grid location
    t =strfind(line,'=');
    
    Riv.IDX(pos,1) = str2double(line(t+1:end));
    % check if IDX is in the list to remove
        line1 = fgetl(fidnml);
        line2 = fgetl(fidnml);
end
fclose (fidnml);
% extract river positions 
for rr=1:length(Riv.IDX)
Riv.lon.data(rr) = Mobj.lon(Riv.IDX(rr,1));
Riv.lat.data(rr) = Mobj.lat(Riv.IDX(rr,1));
Riv.x.data(rr) = Mobj.x(Riv.IDX(rr,1));
Riv.y.data(rr) = Mobj.y(Riv.IDX(rr,1));
end

% %% for each river estimate the annual cycle mean and river
% % flow/concentration relationship
% % rnut = Riv2.N1_p;
% % read tamar nutrient observations. 
% ncfile = fullfile('/data/euryale4/backup/mbe/Experiments/Rosa_rivers/Analysis/tamar_nutrient_nc','gunnislake_dayavg_%s.nc')
% vars = {'si','amm','phs','nitrate','o2'}
% varsR = {'N5_s','N4_n','N1_p','N3_n','O2_o'}
% tamar.time = ncread(sprintf(ncfile,vars{1}),'time');
% for nn=1:length(vars)
% tamar.([varsR{nn},'T']) = ncread(sprintf(ncfile,vars{nn}),'time');
%  tamar.(varsR{nn}) = ncread(sprintf(ncfile,vars{nn}),vars{nn});
% end
%%
Riv.mtime.data = (Riv.time.data+ 678942);
% rr=10;