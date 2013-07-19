function [] = ST_summary(d,varargin)
% Print summary for stats of particle diameter d (m) 
%
% function [] = ST_summary(d,varargin)
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
%      
%
% EXAMPLE USAGE
%    ST_summary(.0005,'temperature',10,'salinity',35,'sdens',2650) 
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% References
%   
%
% Revision history
%   
%==============================================================================

subname = 'ST_summary';  
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

fprintf('       phi          class       d(mm)    Dstar    wset(mm/s)  taucr (Pa)  erate x1e-3(kg/(m^2-s))\n')

phi         = ST_d2phi(d);
phiclass    = ST_wentworth(phi);
Dstar       = ST_Dstar(d);
Wset        = ST_wset(d);
Taucr       = ST_taucr(d);
erate       = ST_erate(d);
fprintf('%10d %20s %8.4f %8.2f %9.4f %8.3f %8.3f\n',phi,phiclass,d*1000,Dstar,Wset*1000,Taucr,1000*erate)


