%
% plot a coastline - quickly with no help as I'm boared of doing this...
%
function [b, arc3] = plot_sms_coast(coast, varargin)

if nargin>1 colour = varargin{1}; else colour = 'g'; end
if nargin>2 plot_type = varargin{2}; else plot_type = 'plot'; end
if nargin>3 I = varargin{3};
else
    I = 1:size(coast.arc,2);
end
if nargin>4 lineW = varargin{4};
else lineW = 0.5;
end

if plot_type(1:4) == 'plot';
    plot_func = @(x, y, col) plot(x, y, 'color', col, 'linewidth', lineW);
elseif plot_type(1:4) == 'patc';
    plot_func = @(x, y, col) patch(x, y, col, 'linewidth', lineW);
elseif plot_type(1:4) == 'm_pl';
    plot_func = @(x, y, col) m_plot(x, y, 'color', col, 'linewidth', lineW);
end

% add in the extra nodes, and use arc2 (modified arc) if it's there...
if isfield(coast, 'arc2')==0
    coast.arc2 = coast.arc;
    coast.node2 = coast.node;
end
for ii=I
    coast.arc3{ii} = [coast.node2{ii}; coast.arc2{ii}; coast.node2{ii}];
end

% plot land masses, one at a time
hold on
for ii=I
    a(ii) = plot_func(coast.arc3{ii}(:,1), coast.arc3{ii}(:,2), colour);
end

b = a(1);
arc3 = coast.arc3;