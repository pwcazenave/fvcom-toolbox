function write_SMS_2dm(file, tri, x, y, bnd)
% Output an SMS 2dm ASCII file from the triangulation given by tri, x and
% y.
%
% write_SMS_2dm(tri, x, y)
% 
% DESCRIPTION:
%   Create an ASCII file in the SMS 2dm format of the triangulation in tri,
%   x and y.
% 
% INPUT:
%   file - file name to save to.
%   tri  - triangulation matrix of the nodes in x and y.
%   x, y - coordinate pairs for the unstructured grid.
%   bnd  - [optional] cell array of open boundary node ids to create node
%          strings in SMS.
% 
% OUTPUT:
%   file - ASCII file in SMS 2dm format.
% 
% EXAMPLE USAGE:
%   write_SMS_2dm('/tmp/test.2dm', Mobj.tri, Mobj.x, Mobj.y)
% 
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
% 
% Revision history:
%   2013-03-11 First version.
% 
%==========================================================================

subname = 'write_SMS_2dm';

global ftbverbose
if ftbverbose
    fprintf('\n'); fprintf(['begin : ' subname '\n']);
end

f = fopen(file, 'w');
if f < 0
    error('Unable to open output file (check permissions?)')
end

if length(x) ~= length(y)
    error('Input coordinate array lengths do not match.')
else
    nn = length(x);
end

nt = size(tri, 1);

% Header
fprintf(f, 'MESH2D\nMESHNAME "mesh"\n');

% Add the connectivity table
for t = 1:nt
    fprintf(f, 'E3T %i %i %i %i 1\n', t, tri(t, :));
end

% Add the list of nodes
for n = 1:nn
    fprintf(f, 'ND %i %.8e %.8e %.8e\n', n, x(n), y(n), 0);
end

% Check we've got some open boundaries and create the relevant node string
% output.
if nargin == 5
    for b = 1:length(bnd)
        c = 0; % counter for the weird nodestring format.

        nodestring = bnd{b}; % current open boundary nodestring indices.
        for ns = 1:length(nodestring)
            % If we're at the end of the nodestring, create a negative node
            % ID. 
            if ns == length(nodestring)
                node_id = -nodestring(ns);
            else
                node_id = nodestring(ns);
            end

            % Increment the counter.
            c = c + 1;
            if c == 1
                % Add the nodestring line prefix and the current node ID.
                fprintf(f, 'NS %i ', node_id);
            elseif c > 0 && c < 10
                fprintf(f, '%i ', node_id);
            elseif c >= 10 || ns == length(nodestring);
                fprintf(f, '%i\n', node_id);
                c = 0;
            end
        end
        % Add a newline at the end of the current nodestring.
        fprintf(f, '\n');
    end
end

% Dump all the (apparently ignored) footer information.
fprintf(f, 'BEGPARAMDEF\nGM  "Mesh"\nSI  0\nDY  0\nTU  ""\nTD  0  0\nNUME  3\nBCPGC  0\nBEDISP  0 0 0 0 1 0 1 0 0 0 0 1\nBEFONT  0 2\nBEDISP  1 0 0 0 1 0 1 0 0 0 0 1\nBEFONT  1 2\nBEDISP  2 0 0 0 1 0 1 0 0 0 0 1\nBEFONT  2 2\nENDPARAMDEF\nBEG2DMBC\nMAT  1 "material 01"\nEND2DMBC\n');

fclose(f);

if ftbverbose
    fprintf('end   : %s\n', subname)
end
