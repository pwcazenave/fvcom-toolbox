% Plot an FVCOM field. This is somewhat similar to plot_field.m but is for
% postprocessing/viewing. It looks for the nv field included in fvcom
% output files. This function runs an animation if the field includes more
% than one time steps.
% 
% INPUT
%   Mobj                    = matlab mesh object 
%   PlotField               = vertex-based field to plot
%   [optional] fid          = the fid of the figure to plot the field in
%   [optional] cli          = the colour limits to use
%   [optional] fig          = make an animated gif - specify filename
%   [optional] axi          = the axis
%
% EXAMPLE USAGE
%    plot_fvcom_field(Mobj,Mobj.zeta,'fid', 1, 'clims', [0 100], fig, 'animation.gif', 'axis', [60000 70000 40000 50000])
%
% Author(s)
%   Rory O'Hara Murray (Marine Scotland Science)
%
function plot_fvcom_field(M, plot_field, varargin)
MJD_datenum = datenum('1858-11-17 00:00:00');

xE = M.x(M.nv)';
yE = M.y(M.nv)';
plot_field = squeeze(plot_field);

if size(plot_field,1)==size(M.nv,1) % plot on elements
    patch_func = @(dummy) patch(xE, yE, dummy', 'linestyle', 'none');
elseif size(plot_field,1)==size(M.x,1) % plot on nodes
    patch_func = @(dummy) patch('Vertices',[M.x, M.y], 'Faces',M.nv, 'Cdata',dummy,'linestyle','none','facecolor','interp');
end

% defaults
clims = [min(plot_field(:)) max(plot_field(:))];
gif = false;

for ii=1:1:length(varargin)
    keyword  = lower(varargin{ii});
    if length(keyword)~=3 continue; end
    switch(keyword(1:3))
        case 'fid' % id of a figure
            fig = varargin{ii+1};
        case 'cli' % colour limits
            clims = varargin{ii+1};
        case 'gif' % make an animated gif
            gif = true;
            gif_filename = varargin{ii+1}
        case 'axi' % axis
            axi = varargin{ii+1};
    end
end

if exist('fig')==0
    fig = figure;
end

for ii=1:size(plot_field,2)
    if ishandle(fig)==0 break; end
    patch_func(plot_field(:,ii));
    colorbar
    set(gca, 'clim', clims)
    title(['time = ' datestr(M.time(ii)+MJD_datenum, 'HH:MM dd/mm/yyyy')])
    
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