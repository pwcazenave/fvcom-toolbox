% example demonstrating reading in an FVCOM mesh and smoothing the bathymetry
%
% function example
%
% DESCRIPTION:
%    Read in FVCOM mesh
%    Read in FVCOM bathymetry
%    Smooth bathymetry
%    Dump new bathymetry
%
% INPUT
%   
% OUTPUT:
%    Smoothed Bathymetry File
%
% EXAMPLE USAGE
%    example2
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================

clear all; close all;

% read the Mesh from an FVCOM grid file
Mobj = read_fvcom_mesh('skg4.3_grd.dat');

% read the bathymetry from an FVCOM grid file
Mobj.h = read_fvcom_bath('skg4.3_dep.dat');  Mobj.have_bath = true;

% smooth bathymetry with 4 iterations of explicit smoother
plot_field(Mobj,Mobj.h,'title','original bathymetry')
Mobj = setup_metrics(Mobj);
[Mobj.h] = smoothfield(Mobj.h,Mobj,0.5,4);
plot_field(Mobj,Mobj.h,'title','smoothed bathymetry');

% dump bathymetry
write_FVCOM_bath(Mobj,'skg4.3_dep_smoothed.dat')

