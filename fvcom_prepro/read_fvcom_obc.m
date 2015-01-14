function Mobj = read_fvcom_obc(Mobj, obcfile) 
% Read fvcom open boundary file 
%
% Mobj = function read_fvcom_obc(obcfile)
%
% DESCRIPTION:
%    Read FVCOM open boundary node specification file.
%
% INPUT:
%   Mobj : existing mesh object into which obc data is added
%   obcfile : fvcom open boundary file
%
% OUTPUT:
%    Mobj.read_obc_nodes : open boundary node cell array (length = number
%    of open boundaries).
%
% EXAMPLE USAGE
%    Mobj = read_fvcom_obc('tst_obc.dat')
%
% Author(s):  
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%    2015-01-12 Create function.
%
%==========================================================================

subname = 'read_fvcom_obc';
global ftbverbose
if ftbverbose
    fprintf('\nbegin : %s\n', subname)
end

have_strings = false;

%--------------------------------------------------------------------------
% read in the FVCOM open boundary nodes
%--------------------------------------------------------------------------
fid = fopen(obcfile,'r');
assert(fid >= 0, 'file: %s does not exist.', obcfile)

reading = true;
while reading
    lin = fgetl(fid);
    if lin ~= -1 % EOF is -1
        lin = strtrim(lin);
        switch lin(1:2)
            case 'OB'
                tmp = regexp(lin, ' = ', 'split');
                nObcNodes = str2double(tmp(end));
                clear tmp
                read_obc_nodes = cell(1, 1);
        end
    else
        reading = false;
    end
end
fclose(fid);

if nObcNodes > 0
    have_strings = true;
    fprintf('nObcNodes = %d\n', nObcNodes)
end

% We have to assume we have only a single open boundary as the _obc.dat
% file has no way of indicating how many open boundaries there are.
fid = fopen(obcfile,'r');
S = textscan(fid, '%f%f%f%*s%*s%[^\n\r]', 'HeaderLines', 1, ...
    'Delimiter', {'\t',' '}, 'MultipleDelimsAsOne', true);
fclose(fid);
read_obc_nodes{1} = S{:, 2};

if have_strings
    Mobj.have_strings = have_strings;
    Mobj.read_obc_nodes = read_obc_nodes;
end

if ftbverbose
    fprintf('\nend   : %s\n', subname)
end