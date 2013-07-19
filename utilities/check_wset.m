function [wset] = ST_wset(d,varargin)
% Calculate settling velocity of particle diameter d (m) in m/s 
%
% function [wset] = ST_wset(d,varargin)
%
% DESCRIPTION:
% Calculate settling velocity of particle diameter d (m) in m/s 
%
% INPUT:
%    d: sediment grain size in m
%    [optional] 'temperature' = temperature of the seawater in C [default=10]
%    [optional] 'salinity'    = salinity of seawater in PSU      [default=35]
%    [optional] 'sdens'       = sediment density in kg/m^3       [default=2650]
%
% OUTPUT:
%    wset: settling velocity in m/s   
%
% EXAMPLE USAGE
%    wset = ST_wset(.0005,'temperature',10,'salinity',35,'sdens',2650) 
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% References
%    Soulsby DMS (102)
%
% Revision history
%   
%==============================================================================

subname = 'ST_wset';  
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


% calculate nu
nu = SW_Kviscosity(T,S);

% calculate rho
dens = SW_Density(T,S);

s = sdens;
dstar = ([grav*(s-1)/(nu^2)])^(1/3)*d;


% calculate wset
wset = (nu/d)*( sqrt(10.36^2 + 1.049*(dstar^3)) - 10.36); 


%fprintf(['end   : ' subname '\n'])
