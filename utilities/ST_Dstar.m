function [dstar] = ST_Dstar(d,varargin)
% Calculate non-dimensional grain size D*
%
% function [dstar] = ST_Dstar(d,varargin)
%
% DESCRIPTION:
%    Convert grain size from d (m) to dimensionless D 
%
% INPUT:
%    d: sediment grain size in m
%    [optional] 'temperature' = temperature of the seawater in C [default=10]
%    [optional] 'salinity'    = salinity of seawater in PSU      [default=35]
%    [optional] 'sdens'       = sediment density in kg/m^3       [default=2650]
%
% OUTPUT:
%    Dstar:  nondimensional grain size
%
% EXAMPLE USAGE
%    dstar = ST_Dstar(.0005,'temperature',10,'salinity',35,'sdens',2650) 
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================

subname = 'ST_Dstar';  
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
                error('incorrect usage of ST_Dstar')
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

% calculate dstar
s = sdens/dens;
dstar = ([grav*(s-1)/(nu^2)])^(1/3)*d;

%fprintf(['end   : ' subname '\n'])
