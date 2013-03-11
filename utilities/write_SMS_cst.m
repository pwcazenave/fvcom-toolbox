function write_SMS_cst(file, x, y)
% Export a set of points to a CST file which can be imported into SMS.
% 
% write_SMS_cst(file, x, y)
% 
% DESCRIPTION:
%   Create an ASCII file in the CST format recognised by SMS from the
%   coordinates specified in (x, y). If x and y are cell arrays, the
%   coordinates in each cell array are treated as different arcs in the CST
%   file (e.g. if you have multiple open boundaries you would like in a
%   single CST file).
% 
% INPUT:
%   file - file name to save to.
%   x, y - coordinate pairs for the open boundary.
% 
% OUTPUT:
%   file - ASCII file in SMS CST format.
% 
% EXAMPLE USAGE:
%   write_SMS_cst(x, y)
% 
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
% 
% Revision history:
%   2013-03-11 First version.
% 
%==========================================================================

subname = 'write_SMS_cst';

global ftbverbose
if ftbverbose
    fprintf('\n'); fprintf(['begin : ' subname '\n']);
end

f = fopen(file, 'w');
if f < 0
    error('Unable to open output file (check permissions?)')
end

if iscell(x) && iscell(y)
    nb = length(x);
else
    nb = 1;
end

% Header
fprintf(f, 'COAST\n');
fprintf(f, '%i\n', nb);

for bb = 1:nb % each boundary
    if iscell(x) && iscell(y)
        np = length(x{bb});
        xx = x{bb};
        yy = y{bb};
    else
        np = length(x);
        xx = x;
        yy = y;
    end

    % The current arc's header
    fprintf(f, '%i\t0.0\n', np);
    
    % All the positions
    for i=1:length(xx);
        fprintf(f, '\t%.6f\t%.6f\t0.0\n', xx(i), yy(i));
    end

end
    
fclose(f);

if ftbverbose
    fprintf('end   : %s \n', subname)
end
