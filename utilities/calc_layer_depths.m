function LayerDepthsFromMSL = calc_layer_depths( M, ncFile, els )
%CALC_LAYER_DEPTHS Find depths of layer centres, including changes over
%time. Output is relative to MSL. 
%FIXME ideally this would be configurable: direction of +ive, and relative
%   to MSL or to seabed.
%   Inputs: M: mesh object. Must contain x,y,tri,h (populated).
%           ncFile: Char or string. Name of FVCOM netCDF output file containing siglay_center field
%           els: vector. Element numbers of interest.
%   
%   Output: LayerDepths: Double matrix of dimensions layer x timestep x
%                           element
%                        Gives the depth (+ive down from MSL) at the
%                        centre of each layer at the centre of each element.

% Simon Waldman / PNNL, May 2019.

global ftbverbose;
if ftbverbose
    [~, subname] = fileparts(mfilename('fullpath'));
    fprintf('\nbegin : %s\n', subname)
end

%check inputs.
assert( isstruct(M) && all( isfield( M, {'x', 'y', 'tri', 'h'} ) ), ...
    'Mobj must be a struct containing x,y,tri and h fields' );
assert( M.have_xy, 'M.have_xy is false.' ); 
ncFile = convertStringsToChars(ncFile);
assert( exist(ncFile, 'file') == 2, 'Can''t find ncFile %s.', ncFile );
assert( max( M.h ) > 0, 'M.h appears to contain all zeroes. Possibly you need to import bathymetry.' );
assert( isvector(els), 'els should be a vector.');
NumEls = length(els);

% load the layer distribution. This gives the depth of each layer in each
% element as a proportion of the total depth at that location. Dims are
% element x layer.
siglayc = ncread(ncFile, 'siglay_center');
NumLayers = size(siglayc, 2);

%if the M object already has a hc field, use it. If not, calculate it.
if ~isfield( M, 'hc' ) || max( M.hc ) == 0
    if ftbverbose
        disp('Calculating hc from h. NB If x and y are actually lon and lat, this will be inaccurate.');
    end
    M.hc = mean( M.h( M.tri ),2 );
end
assert( all(els <= M.nElems), 'One or more requested element numbers is higher than the number of elements in the mesh.');

if ftbverbose
    disp('Loading free surface elevations from ncFile...');
end
% load the free surface changes from the ncfile.
zeta = ncread( ncFile, 'zeta' ); %this has dims of node x timestep
if ftbverbose
    disp('done.');
end
NumTS = size( zeta, 2 );

%convert to elements (that we need) & calc total & layer depths
LayerDepthsFromMSL = nan( NumLayers, NumTS, NumEls );
for e = 1:NumEls
    el = els(e);
    el_zeta = mean( zeta(M.tri(el,:),:), 1 );
    TotalDepths = el_zeta + repmat( M.hc(el), 1, NumTS ); %1xNumTS vector of depths
    LayerDepths = repmat( TotalDepths, NumLayers, 1 ) .* repmat( siglayc(el,:)', 1, NumTS ) * -1;
    % the -1 is needed because siglayc is negative.
    LayerDepthsFromMSL(:,:,e) = LayerDepths - repmat( el_zeta, 10, 1 );
end

end

