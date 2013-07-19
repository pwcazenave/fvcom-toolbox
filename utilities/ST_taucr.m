function [taucr] = ST_taucr(d,varargin)
% Calculate critical shear stress in Pascals 
%
% function [wset] = ST_taucr(d,varargin)
%
% DESCRIPTION:
% Calculate critical shear stress for threshold of motion in Pa
%
% INPUT:
%    d: sediment grain size in m
%    [optional] 'temperature' = temperature of the seawater in C [default=10]
%    [optional] 'salinity'    = salinity of seawater in PSU      [default=35]
%    [optional] 'sdens'       = sediment density in kg/m^3       [default=2650]
%
% OUTPUT:
%    taucr:  critical shear stress in N/m^2  
%
% EXAMPLE USAGE
%    TCR = ST_taucr(.0005,'temperature',10,'salinity',35,'sdens',2650) 
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% References
%    Soulsby Dynamics of Marine Sands (SC77)
%
% Revision history
%   
%==============================================================================

subname = 'ST_taucr';  
%fprintf('\n')
%fprintf(['begin : ' subname '\n'])

% constants 
grav  = 9.8106;   %g
T     = 10;       %T (C)
S     = 35;       %S (PSU)
sdens = 2650;     %sediment density in kg/m^3

% parse arguments
for i=1:2:length(varargin)-1
        keyword  = lower(varargin{i});
        if( ~ischar(keyword) )
                error('incorrect usage of ST_wset')
        end;

        switch(keyword(1:3))

        case 'tem'
             T = varargin{i+1};
        case 'sal'
             S = varargin{i+1};
        case 'sde'
             sdens = varargin{i+1}; 
        otherwise
                error(['Can''t understand value for:' keyword]);
        end; %switch keyword
end;


% calculate rho
dens = SW_Density(T,S);

% calculate dstar
dstar = ST_Dstar(d,'temp',T,'sal',S,'sdens',sdens);

% calculate theta_cr
theta_cr = (0.30/(1+1.2*dstar)) + 0.055*[1 - exp(-.020*dstar)];

% calculate taucr
taucr = theta_cr*grav*(sdens-dens)*d;



%fprintf(['end   : ' subname '\n'])
