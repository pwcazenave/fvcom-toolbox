function  [X,Y,BUTTON,SCALEMAT] = ginput2(varargin)
%GINPUT2   Graphical input from mouse with zoom, pan, plot and scaling.
%   
%   SYNTAX:
%                        XY = ginput2;
%                        XY = ginput2(DoScale);          true or false
%                        XY = ginput2(...,PlotOpt);      '.r' for example
%                        XY = ginput2(...,'KeepZoom');   vs. 'UnZoom'
%                        XY = ginput2(N,...);
%                        XY = ginput2(...);
%                     [X,Y] = ginput2(...);
%              [X,Y,BUTTON] = ginput2(...);
%     [X,Y,BUTTON,SCALEMAT] = ginput2(...);
%
%   INPUT:
%     DoScale    - Single logical specifying whether the IMAGE should be
%                  interactively scaled (georeferenced), or it can be the
%                  2x4 SCALEMAT matrix for automatically scaling.
%                  DEFAULT: false (do not scales/georeferences)
%     PlotOpt    - String and/or parameter/value pairs specifying the drawn
%                  points optional inputs (see PLOT for details). 
%                  DEFAULT: 'none' (do not plot any point)
%     'KeepZoom' - When finishing selection by default the zoom is
%                  restored. By using this option this feature is ignored.
%                  DEFAULT: 'UnZoom' (restores original axis limits)
%     N          - Number of points to be selected. One of 0,1,2,...,Inf
%                  DEFAULT: Inf (selects until ENTER or ESCAPE is pressed)
%
%   OUTPUT:
%     XY        - [X(:) Y(:)] axis coordinate(s).
%     X         - X-coordinate(s).
%     Y         - Y-coordinate(s).
%     BUTTON    - Last pressed button.
%     SCALEMAT  - 2x4 matrix specifying the coordinates of two different
%                 points (1) and (2) in the Image coordinates (pixels) and
%                 the User coordinates (data):
%                                                Point 1     Point 2
%                   Image coord (pixels):     [ (I1x,I1y)   (I2x,I2y) 
%                   User  coord (data)  :       (U1x,U1y)   (U2x,U2y) ]
%                 to be use for scaling/georeferencing.
%
%   DESCRIPTION:
%     This program uses MATLAB's GINPUT function to get the coordinates
%     of a mouse-selected point in the current figure (see GINPUT for
%     details), but with five major improvements:
%                  1. ZOOMING  (left click)
%                  2. PANNING  (dragging mouse)
%                  3. DELETING (last selected point)
%                  4. PLOTING  (temporarily the selected points)
%                  5. SCALING or GEOREFERENCE IMAGES.
%     The differences are:
%      a) Obviously, the SCALEOPT and PlotOpt optional arguments.
%      b) When click is made outside the axes, it is ignored.
%      c) When LEFT-click, ZOOM-IN is performed right into the selected
%         point (PANNING).
%      d) When RIGHT-click, the point is selected (normal).
%      e) When DOUBLE-click, ZOOM-OUT is done.
%      f) When MIDDLE-click, ZOOM-RESET is done (see ZOOM for details).
%      g) When dragging while pressed left-click PAN is done (until the
%         button is released).
%      h) When pressed any KEY follows the next rules: 
%          A) If ENTER is pressed, the selection is terminated. If no point
%             was already selected, the outputs are empty's.
%          B) If BACKSPACE key is pressed, the last selected point is
%             deleted and the selection continues.
%          C) If SPACEBAR the mouse current position or NANs coordinates
%             are saved, depending whether the mouse was inside or outside
%             any of the current figure axes, respectively. In this latter
%             case, the selection is NOT counted as one of the N points.
%             Besides, when drawing the color is changed. Then, the outputs
%             may not be of length N.
%
%   NOTE:
%     * Optional inputs use its DEFAULT value when not given or [].
%     * Optional outputs may or not be called.
%     * String inputs may be shortened, as long as they are unambiguous.
%       Case is ignored.
%     * The function can be used for interactively digitalize/vectorize
%       RASTER images with:
%       >> ginput(true)
%     * The function can be used only as a georeference function with 
%       >> ginput2(0,true)
%     * The scale/georeference only works when the current axes has an
%       IMAGE type children (see Image for details). 
%     * The x and y data from axes and image are changed when scale/
%       georeference is used.
%     * The drawn points are deleted from the graphics once the selection
%       is finished. 
%     * The priority of the inputs are: N, then SCALEOPT and finally
%       PlotOpt. If the first (integer) is missing, the next is taken into
%       account (logical or 2x4 matrix) and so on.
%
%   EXAMPLE:
%     % Selects until ENTER is pressed:
%         xy = ginput2;
%     % Selects 5 points:
%         [x,y] = ginput2(5);
%     % Gets pressed button:
%         [x,y,button] = ginput2(1);
%     % Scales image and select 4 points temporarily coloring them in
%       black. Besides to not ZOOM OUT at the end:
%         imagesc(peaks(40))
%         [x,y,button,scalemat] = ginput2(4,true,'k*','KeepZoom');
%         hold on, plot(x,y,'or'), hold off
%
%   SEE ALSO:
%     GINPUT, PLOT.
%
%
%   ---
%   MFILE:   ginput2.m
%   VERSION: 3.1 (Nov 12, 2009) (<a href="matlab:web('http://www.mathworks.com/matlabcentral/fileexchange/authors/11258')">download</a>) 
%   MATLAB:  7.7.0.471 (R2008b)
%   AUTHOR:  Carlos Adrian Vargas Aguilera (MEXICO)
%   CONTACT: nubeobscura@hotmail.com
