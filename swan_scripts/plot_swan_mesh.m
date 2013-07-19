function plot_swan_mesh(bathfile,nodefile,gridfile)
% Plot an unstructured SWAN mesh
%
% function plot_swan_mesh(bathfile,nodefile,gridfile)
%
% DESCRIPTION:
%    plot an unstructure SWAN mesh and bathymetry
%
% INPUT 
%   bathfile = SWAN bathymetry file
%   nodefile = SWAN vertex file
%   gridfile = SWAN connectivity file    
%
% OUTPUT:
%    plot of Mesh and bathymetry
%
% EXAMPLE USAGE
%    plot_swan_mesh('tst.bot','tst.node','tst.ele') 
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================

subname = 'plot_swan_mesh';
fprintf('\n')
fprintf(['begin : ' subname '\n'])

%------------------------------------------------------------------------------
% read in the SWAN bathymetry, connectivity, node files
%------------------------------------------------------------------------------


[i,i1,i2,i3] = textread(gridfile,'%d %d %d %d','headerlines',1);
tri(:,1) = i1;
tri(:,2) = i2;
tri(:,3) = i3;

[i,x,y,type] = textread(nodefile,'%d %f %f %d\n','headerlines',1);
x(:,1) = x;
x(:,2) = y;

[bath] = textread(bathfile,'%f\n');

patch('Vertices',x,'Faces',tri,...
       'Cdata',bath,'edgecolor','interp','facecolor','interp');
colorbar

hold on
bnodes = find(type==1);
plot(x(bnodes,1),y(bnodes,1),'k+')
obnodes = find(type==2);
plot(x(obnodes,1),y(obnodes,1),'ro')


fprintf(['end   : ' subname '\n'])

 
