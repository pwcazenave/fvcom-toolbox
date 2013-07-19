function Y = runmean(X, m, dim, modestr) ;
% RUNMEAN - Very fast running mean (aka moving average) filter
%   For vectors, Y = RUNMEAN(X,M) computes a running mean (also known as
%   moving average) on the elements of the vector X. It uses a window of
%   2*M+1 datapoints. M an positive integer defining (half) the size of the
%   window. In pseudo code: 
%     Y(i) = sum(X(j)) / (2*M+1), for j = (i-M):(i+M), and i=1:length(X) 
%
%   For matrices, Y = RUNMEAN(X,M) or RUNMEAN(X,M,[]) operates on the first   
%   non-singleton dimension of X. RUNMEAN(X,M,DIM) computes the running
%   mean along the dimension DIM.
%
%   If the total window size (2*M+1) is larger than the size in dimension
%   DIM, the overall average along dimension DIM is computed.
%
%   As always with filtering, the values of Y can be inaccurate at the
%   edges. RUNMEAN(..., MODESTR) determines how the edges are treated. MODESTR can be
%   one of the following strings:
%     'edge'    : X is padded with first and last values along dimension
%                 DIM (default)
%     'zero'    : X is padded with zeros
%     'mean'    : X is padded with the mean along dimension DIM 
%
%   X should not contains NaNs, yielding an all NaN result. NaNs can be
%   replaced by using, e.g., "inpaint_nans" created by John D'Errico.
%
%   Examples
%     runmean([1:5],1) 
%       % ->  1.33  2  3  4 4.67
%     runmean([1:5],1,'mean') 
%       % ->  2 2 3 4 4
%     runmean([2:2:10],1,1) % dimension 1 is larger than 2*(M=1)+1 ...
%       % -> 2 4 6 8 10
%     runmean(ones(10,7),3,2,'zero') ; % along columns, using mode 'zero'
%     runmean(repmat([1 2 4 8 NaN 5 6],5,1),2,2) ; 
%       % -> all NaN result
%     A = rand(10,10) ; A(2,7) = NaN ;
%     runmean(A,3,2) ; 
%       % -> column 7 is all NaN
%     runmean(1:2:10,100) % mean
%       % -> 5 5 5 5 5
%
%   This is an incredibly fast implementation of a running mean, since
%   execution time does not depend on the size of the window.
%
%   See also MEAN, FILTER

% for Matlab R13
% version 3.0 (sep 2006)
% Jos van der Geest
% email: jos@jasen.nl

% History:
%   1.0 (2003) created, after a snippet from Peter Acklam (?)
%   1.1 (feb 2006) made suitable for the File Exchange (extended help and
%       documentation)
%   1.2 (feb 2006) added a warning when the window size is too big
%   1.3 (feb 2006) improved help section
%   2.0 (sep 2006) working across a dimension of a matrix. 
%   3.0 (sep 2006) several treatments of the edges. 

% Acknowledgements: (sep 2006) Thanks to Markus Hahn for the idea of
% working in multi-dimensions and the way to treat edges. 

error(nargchk(2,4,nargin)) ;

if ~isnumeric(m) || (numel(m) ~= 1) || (m < 0) || fix(m) ~= m,
    error('The window size (M) should be a positive integer') ;
end

if nargin == 2,
    dim = [] ;
    modestr = 'edge' ;
elseif nargin==3,
    if ischar(dim),
        % no dimension given
        modestr = dim ;
        dim = [] ;
    else 
        modestr = 'edge' ;
    end
end

modestr = lower(modestr) ;

% check mode specifier
if ~ismember(modestr,{'edge','zero','mean'}),
    error('Unknown mode') ;
end

szX = size(X) ;
if isempty(dim),
    dim = min(find(szX>1)) ;
end

if m == 0 || dim > ndims(X),
    % easy
    Y = X ;
else
    mm = 2*m+1 ;
    if mm >= szX(dim),
        % if the window is larger than X, average all
        sz2 = ones(size(szX)) ; 
        sz2(dim) = szX(dim) ;
        Y = repmat(mean(X,dim),sz2) ;
    else
        % here starts the real stuff
        % shift dimensions so that the desired dimensions comes first
        [X, nshifts] = shiftdim(X, dim-1); 
        szX = size(X) ;
        % make the rest of the dimensions columns, so we have a 2D matrix
        % (suggested of Markus Hahn)
        X = reshape(X,szX(1),[]) ; 
        % select how to pad the matrix
        switch (modestr),
            case 'edge'
                % pad with first and last elements
                Xfirst = repmat(X(1,:),m,1) ;
                Xlast = repmat(X(end,:),m,1) ;
            case 'zero'
                % pad with zeros
                Xfirst = zeros(m,1) ;
                Xlast= zeros(m,1) ;
            case 'mean',
                % pad with the average
                Xfirst = repmat(mean(X,1),m,1) ;
                Xlast = Xfirst ;
        end        
        % pad the array
        Y = [zeros(1,size(X,2)) ; Xfirst ; X ; Xlast] ;       
        % the cumsum trick (by Peter Acklam ?)
        Y = cumsum(Y,1) ;
        Y = (Y(mm+1:end,:)-Y(1:end-mm,:)) ./ mm ;
        
        % reshape into original size
        Y = reshape(Y,szX)   ;
        % and re-shift the dimensions
        Y = shiftdim(Y,ndims(Y)-nshifts) ;
    end
end

% =====================
%  CODE OF VERSION 1.3 
% =====================

% function Y = runmean(X,m) ;
% % RUNMEAN - Very fast running mean filter for vectors
% %   Y = RUNMEAN(X,M) computes a running mean on vector X using a window of
% %   2*M+1 datapoints. X is a vector, and M an positive integer defining
% %   (half) the size of the window. In pseudo code:
% %     Y(i) = sum(X(j)) / (2*M+1), for j = (i-M):(i+M), and i=1:length(X)
% %
% %   If the total window size (2M+1) is larger than the length of the vector, the overall
% %   average is returned.
% %
% %   Example:
% %     runmean(1:10,1) % ->
% %     [1.3333 2 3 4 5 6 7 8 9 9.6667]
% %
% %   This is an incredibly fast implementation of a running average, since
% %   execution time does not depend on the size of the window.
% %
% %   X should not contains NaNs (a NaN will result in a all NaN result)
% %   At both ends the values of Y can be inaccurate, as the first and last
% %   values of X are used multiple times. 
% %
% %   See also MEAN
% 
% % for Matlab R13
% % version 1.3 (feb 2006)
% % Jos van der Geest
% % email: jos@jasen.nl
% 
% % History:
% % 1.0 (2003) created, after a snippet from Peter Acklam (?)
% % 1.1 (feb 2006) made suitable for the File Exchange (extended help and
% % documentation)
% % 1.2 (feb 2006) added a warning when the window size is too big
% % 1.3 (feb 2006) improved help section
% 
% error(nargchk(2,2,nargin)) ;
% 
% sz = size(X) ;
% 
% if numel(sz) ~= 2 || (min(sz) ~= 1),
%     error('X should be a vector') ;
% end
% 
% if any(isnan(X)),
%     error('NaNs cannot be dealt with') ;
% end
% 
% if ~isnumeric(m) || (numel(m) ~= 1) || (m < 0) || fix(m) ~= m,
%     error('The window size (M) should be a positive integer') ;
% elseif m == 0,
%     Y = X ;
%     return ;
% end
% 
% mm = 2*m+1 ;
% 
% if mm >= prod(sz),
%     % if the window is larger than X, average all
%     warning('Window size is larger than the length of the vector.')
%     Y = repmat(mean(X),sz) ;
% else
%     % the cumsum trick ...
%     Y = [repmat(X(1),m,1) ; X(:) ; repmat(X(end),m,1)] ;
%     Y = [0 ; cumsum(Y)] ;
%     Y = (Y(mm+1:end)-Y(1:end-mm)) / mm ;
%     Y = reshape(Y,sz) ;
% end

