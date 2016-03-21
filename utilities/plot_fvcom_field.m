% Plot an FVCOM field. This is somewhat similar to plot_field.m but is for
% postprocessing/viewing. It looks for the nv field included in fvcom
% output files. This function runs an animation if the field includes more
% than one time steps.
%
% plot_fvcom_field(Mobj, PlotField, 'fid', figure_id, 'cli', colour_lims, 'gif',
% filename, 'axi', axis_range, 'pll', 'grd', colour);
%
% INPUT
%   Mobj                    = matlab mesh object with the following fields:
%       - lon, lat, x, y    = nodal positions (spherical and/or cartesian)
%       - nv, tri           = connectivity table (called either nv or tri)
%       - [optional] time   = Modified Julian Day time series
%   PlotField               = vertex-based field to plot
%   [optional] fid          = the fid of the figure to plot the field in - specify figure id
%   [optional] cli          = the colour limits to use - specify the limits
%   [optional] gif          = make an animated gif - specify filename
%   [optional] axi          = the axis - specify axis range
%   [optional] pll          = the axis
%   [optional] grd          = add gridlines - specify colour
%
% EXAMPLE USAGE
%    plot_fvcom_field(Mobj, Mobj.zeta, 'fid', 1, 'cli', [0 100], 'gif', 'animation.gif', 'axi', [60000 70000 40000 50000])
%
% Author(s)
%   Rory O'Hara Murray (Marine Scotland Science)
%
% Developments:
% 2014-05-22: Changed the way fig id is checked, not using 'exist' anymore.
% 2014-08-15: Added the axis command in
%
function [a] = plot_fvcom_field(M, plot_field, varargin)
MJD_datenum = datenum('1858-11-17 00:00:00');

% check to see if nv or tri should be used.
if isfield(M, 'nv')
    nv = M.nv;
elseif isfield(M, 'tri')
    nv = M.tri;
end

% check to see if a time variable is there or not
if isfield(M, 'time') %& size(M.time,1)>1
    time_flag = true;
else
    time_flag = false;
end

% defaults
clims = [min(plot_field(:)) max(plot_field(:))];
if clims(1)==clims(2)
    clims(1)=clims(1)-0.1;
    clims(2)=clims(2)+0.1;
end
gif = false;
grd = false;
plot_ll = false;
fig_flag = false;
axis_flag = false;
title_flag = false;
legend_text_flag = false;
quiver_flag = false;
quiver2_flag = false;

for ii=1:1:length(varargin)
    keyword  = lower(varargin{ii});
    if length(keyword)~=3
        continue
    end
    switch keyword
        case 'fid' % id of a figure
            fig = varargin{ii+1};
            fig_flag = true;
        case 'cli' % colour limits
            clims = varargin{ii+1};
        case 'gif' % make an animated gif
            gif = true;
            gif_filename = varargin{ii+1};
        case 'axi' % axis
            axis_flag = true;
            axi = varargin{ii+1};
        case 'grd' % grid lines
            grd = true;
            edgecolor = varargin{ii+1};
        case 'pll'
            plot_ll = true;
        case 'tit'
            title_flag = true;
            fig_title = varargin{ii+1};
        case 'leg'
            legend_text_flag = true;
            legend_text = varargin{ii+1};
        case 'qui'
            quiver_flag = true;
            quiverData = varargin{ii+1};
        case 'qu2'
            quiver2_flag = true;
            quiverData = varargin{ii+1};
    end
end

if plot_ll
    x = M.lon;
    y = M.lat;
else
    x = M.x;
    y = M.y;
end

if not(axis_flag)
    axi = [min(x) max(x) min(y) max(y)];
end

xE = x(nv)';
yE = y(nv)';
plot_field = squeeze(plot_field);

if size(plot_field,1)==size(nv,1) % plot on elements
    if grd
        patch_func = @(dummy) patch(xE, yE, dummy', 'edgecolor', edgecolor);
    else
        patch_func = @(dummy) patch(xE, yE, dummy', 'linestyle', 'none');
    end
elseif size(plot_field,1)==size(x,1) % plot on nodes
    if grd
        patch_func = @(dummy) patch('Vertices',[x, y], 'Faces',nv, 'Cdata',dummy,'edgecolor', edgecolor,'facecolor','interp');
    else
        patch_func = @(dummy) patch('Vertices',[x, y], 'Faces',nv, 'Cdata',dummy,'linestyle','none','facecolor','interp');
    end
end

if not(fig_flag)
    fig = figure;
end

for ii=1:size(plot_field,2)
    if ishandle(fig)==0 break; end
    a = patch_func(plot_field(:,ii));
    c = colorbar;
    if legend_text_flag set(get(c, 'ylabel'), 'string', legend_text); end
    set(gca, 'clim', clims);
    axis(axi)
    if title_flag
        title(fig_title)
    elseif time_flag
        title(['time = ' datestr(double(M.time(ii))+MJD_datenum, 'HH:MM dd/mm/yyyy')])
    end
    if quiver_flag
        hold on
        quiver(quiverData.X, quiverData.Y, quiverData.U(:,:,ii), quiverData.V(:,:,ii), 'w');
        hold off
    elseif quiver2_flag
        hold on
        quiver(quiverData.X, quiverData.Y, quiverData.U(:,ii), quiverData.V(:,ii), 'color', 'w')
        hold off
    end

    if gif
        axis off
        set(gcf, 'color', 'w')
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if ii == 1;
            imwrite(imind,cm,gif_filename,'gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,gif_filename,'gif','WriteMode','append');
        end
    else
        pause(0.01);
    end
end

return
