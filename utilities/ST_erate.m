function [erate] = ST_erate(d,varargin)
% Calculate erosion rate in kg/(m^2-s)
%
% function [erate] = ST_erate(d,varargin)
%
% DESCRIPTION:
% Calculate erosion rate of sediment diameter d (m) in kg/(m^2-s)
%   See NOTE below
%
% INPUT:
%    d: sediment grain size in m
%    [optional] 'temperature' = temperature of the seawater in C [default=10]
%    [optional] 'salinity'    = salinity of seawater in PSU      [default=35]
%    [optional] 'sdens'       = sediment density in kg/m^3       [default=2650]
%
% OUTPUT:
%    erate in kg/(m^2-s)
%
% EXAMPLE USAGE
%    erate = ST_erate(.0005,'temperature',10,'salinity',35,'sdens',2650) 
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% References
%    Blaas etal, Cont. Shelf. Res., 27, 2007
%
% Note
%    THIS CALCULATION IS A HACK BASED ON A CURVE FIT OF SEVERAL STUDIES
%    Need to modify to use formula of Drake and Cacchione, CSR 9, 1989
%    or other.  Someone fix this if possible.
%
% Revision history
%   
%==============================================================================

subname = 'ST_erate';  
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

% calculate settling velocity
wset = ST_wset(d,'temperature',T,'salinity',S,'sdens',sdens);

% calculate erosion rate
erate = 2.666e-4*wset*1000. - 2.51e-9*sdens;



%fprintf(['end   : ' subname '\n'])
