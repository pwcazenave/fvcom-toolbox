function data = read_sms_map(file)
% Reads a .map file from SMS
%
%  data = read_sms_map(filename)
%
% DESCRIPTION:
%   Reads a map file from SMS consisting of node coordinates for boundary
%   lines. Each different arc is read into a different cell array (e.g. if
%   you have multiple islands). Ensure that the polygons are built in SMS 
%   to make handling islands easier.
%
% INPUT:
%   file - file name to read from.
%
% OUTPUT:
%   data.arc - cell array of coordinate pairs for the boundaries and
%       elevation
%   data.arcID - cell array of SMS ARC IDs
%   data.arcN - cell array of number of nodes in each arc cell array
%
% EXAMPLE USAGE:
%  arc = read_sms_map('/home/user/data/sms_project_v1.map')
%
% Author(s):
%   Ricardo Torres (Plymouth Marine Laboratory)
%
% Revision history:
%   2013-05-03 First version.
%
%==========================================================================

subname = 'read_sms_map';

global ftbverbose
if ftbverbose
    fprintf('\n'); fprintf(['begin : ' subname '\n']);
end
f = fopen(file, 'rt');
if f < 0
    error('Unable to open output file (check permissions?)')
end
frewind(f)
loopnode = 0;
looparc = 0;
data = [];
while ~feof(f)
    % Read NODES found
    line=fgetl(f);
    while ~(strcmpi(line,'ARC') | strcmpi(line,'NODE'))
        
        line=fgetl(f);
        if  feof(f);fclose (f); return; end
    end
    switch line
        case 'NODE'
            
            loopnode = loopnode + 1;
            line = fgetl(f);
            data.node(loopnode) = textscan(line, '%*s %f %f %*f', 'CollectOutput', 1);
            line = fgetl(f);
            data.nodeID(loopnode) = textscan(line, '%*s %u');
            line = fgetl(f);
            
        case 'ARC'

            looparc = looparc + 1;
            
            line=fgetl(f);
            data.arcID{looparc} = textscan(line, '%*s %u');
            % skip two lines
            dump = fgetl(f); line = fgetl(f);
            data.arcnode(looparc) = textscan(line, '%*s %u %u', 'CollectOutput', 1);
            line = fgetl(f);
            if strcmpi(line, 'END') ; break; end
            % read number of vertices in ARC
            data.arcN(looparc) = textscan(line, '%*s %u');
            data.arc(looparc) = textscan(f, '%f %f %*f', data.arcN{looparc}, 'CollectOutput', 1, 'Delimiter', '\n');
        otherwise
            
    end
end
fclose(f)

if ftbverbose;
    fprintf(['end   : ' subname '\n'])
end

return
