function [Mobj] = get_POLPRED_spectide(Mobj, POLPRED)
% Extract tidal harmonic phases and amplitudes from POLPRED ASCII files.
%
% get_POLPRED_spectide(Mobj, POLPRED)
%
% DESCRIPTION:
%    Extract POLPRED harmonic amplitude and phases for the nearest point in
%    the POLPRED grid to the open boundary nodes in Mobj.
%
% INPUT:
%   Mobj    = MATLAB mesh object (see read_sms_mesh.m)
%   POLPRED = ASCII file of the POLPRED harmonics
%
% OUTPUT:
%    Mobj  = MATLAB mesh object with two new arrays:
%       phase - Harmonic phases at each open boundary point
%       amp   - Harmonic amplitudes at each open boundary point
%
% EXAMPLE USAGE
%    Mobj = get_POLPRED_spectide(Mobj, '/path/to/POLPRED.txt')
%
% Author(s):  
%    Pierre Cazenave (Plymouth Marine Laboratory)
% 
% Revision history
%    2012-11-15 First version. Based in part on tide_tools.py from the
%    fvcom-py Python toolbox (https://bitbucket.org/pwcazenave/fvcom-py)
%    and Ricardo Torres' searchtides.m.
%   
%==========================================================================

subname = 'get_POLPRED_spectide';

global ftbverbose;
if ftbverbose
  fprintf('\n')
  fprintf(['begin : ' subname '\n'])
end

% Check we have spherical coordinates in Mobj (we can't extract harmonics
% without them (easily)).
if ~Mobj.have_lonlat
    error('Spherical coordinates absent from Mobj. Cannot extract harmonics from cartesian coordinates.')
end

% Read the POLPRED header into a struct of header names plus their values.
fid = fopen(POLPRED,'rt');
if(fid < 0)
	error(['file: ' POLPRED ' does not exist']);
end

if ftbverbose
  fprintf(['reading from: ' POLPRED '\n'])
  fprintf('extracting header\n')
end

readingHeader = true;
header = struct();
while readingHeader
	lin = fgetl(fid);
    if isempty(lin)
        % Got to the end of the header
        readingHeader = false;
    else
        % We have some header information. Load it into a struct.
        key = regexp(lin, ':', 'split');
        header.(strtrim(regexprep(key{1}, ' ', '_'))) = strtrim(key{2});
    end
end

% Reformat the list of harmonics into a more sensible format
header.Harmonics = regexp(header.Harmonics, ' ', 'split');

% Get the positions in header.Harmonics for the constituents in which we're
% interested.

pInd = 1:length(header.Harmonics);
pIndUse = nan(length(Mobj.Components), 2);
for i=1:length(Mobj.Components)
    pPos = pInd(strcmp(Mobj.Components{i}, header.Harmonics));
    if isempty(pPos)
        warning('Supplied constituent (%s) is not present in the POLPRED data', Mobj.Components{i}) %#ok<WNTAG>
    else
        % Make index start at zero so the multiplication works, but
        % compensate for that once the offset has been applied. Also add
        % offset for the 2 columns (amplitude and phase).
        pIndUse(i, :) = (repmat((pPos - 1 ) * 6, 1, 2) + 1) + (0:1);
    end
end
% Add three to offset by the lat, lon and flag columns
pIndUse = pIndUse + 3;

% Now we're at the data. Load it all into a massive array.
if ftbverbose
  fprintf('extracting data\n')
end

readingData = true;
i = 0;
% Preallocate data to something big and then cut back afterwards (if
% necessary). Get the number of columns from the header and multiply by 6
% (amplitude and phases for z, u and v). Add three for the lat, lon and
% flag columns). The rows is the number of data lines in my files for the
% northwest European shelf domain.
nCols = 3 + (str2double(header.Number_of_harmonics) * 6);
data = nan(358797, nCols);
if ftbverbose
    tic
end
while readingData
    lin = fgetl(fid);
    if lin ~= -1 % EOF is -1
        i = i + 1;
        if ftbverbose
            if mod(i, 10000) == 0
                fprintf('line %i\n', i)
            end
        end
        % str2double doesn't work without a couple of calls to regexp,
        % which makes it ~20x slower than str2num on its own. The regexp
        % approach is still here if you don't believe me.
        data(i, :) = str2num(strtrim(lin)); %#ok<ST2NM>
%         data(i, :) = str2double(regexp(regexprep(strtrim(lin), '\s+', ' '), ' ', 'split'));
    else
        if ftbverbose
            fprintf('end of file at line %i\n', i)
        end
        readingData = false;
    end
end
if ftbverbose
    toc
end

fclose(fid);

% Clear out any additional NaNs in data from preallocation.
data = reshape(data(~isnan(data)), i, nCols);

% Now we have the data, identify the indices of data which correspond to
% the nearest point to each open boundary point. This approach may not be
% the best: it might instead be better to simply read in the positions and
% create an index which we use to extract the harmonics of interest.
% However, we've got this far so might as well carry on.

% Get the open boundary node IDs with which to extract their locations
tmpObcNodes = Mobj.obc_nodes';
ObcNodes = tmpObcNodes(tmpObcNodes~=0)';
obc_lon = Mobj.lon(ObcNodes);
obc_lat = Mobj.lat(ObcNodes);

% For each position, find the nearest POLPRED value. Use the
% find_nearest_pt.m logic to get the nearest point (we can't use the
% function here because the values for which we're searching aren't in
% Mobj).
distance = nan(size(obc_lon));
point = nan(size(distance));
% Omit the NaNs in the indices from POLPRED when calculating the output
% array size.
amp = nan(length(obc_lon), length(pIndUse(~isnan(pIndUse(:, 1)), 1)));
phase = nan(size(amp));
for i=1:length(obc_lon)
    radvec = sqrt((obc_lon(i)-data(:,2)).^2 + (obc_lat(i)-data(:,1)).^2);
    [distance(i), point(i)] = min(radvec);
    % Get the amplitude and phase for each constituent (in order of
    % Mobj.Components). Check for and omit NaNs here (for missing tidal
    % constituents in the supplied list and what's given in POLPRED).
    amp(i, :) = data(point(i), pIndUse(~isnan(pIndUse(:, 1)), 1));
    phase(i, :) = data(point(i), pIndUse(~isnan(pIndUse(:, 1)), 2));
end

% Check for and warn about NaNs (-999.9 in POLPRED data).
if sum(amp(:)==-999.9) > 0
    warning('NaN values (-999.9 in POLPRED terms) in the amplitude data. Are your boundaries on land?') %#ok<WNTAG>
end
if sum(phase(:)==-999.9) > 0
    warning('NaN values (-999.9 in POLPRED terms) in the phase data. Are your boundaries on land?') %#ok<WNTAG>
end

Mobj.amp = amp;
Mobj.phase = phase;

% Plot the open boundary positions and the closest POLPRED point.
% figure(1000)
% plot(obc_lon, obc_lat, 'o')
% hold on
% plot(data(point,2), data(point,1), 'rx')

