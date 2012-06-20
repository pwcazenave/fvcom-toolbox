function [tauw] = calc_tauw(z0,ncfile,writenetcdf) 

% DESCRIPTION:
%    Function used to calculate Wave Bed Shear Stress, tauw, from the
%    bottom orbital velocity (Ubot) and the bottom wave period (TmBot). 
%    This function is currently set up to extract these values from SWAN
%    output. Currently temperature, salinity and z0 are set to be constant
%    accross the domain. Using actual data could improve accuaracy.
%    
%    REQUIRES FVCOM_TOOLBOX to be in your MATLAB path
% 
%
% INPUT 
%   z0           = bed roughness length (z0=d50/12) [mm]
%   <=={possibly use constant value 0.05 Soulsby pg. 49}==>
%   
%   ncfile       = netcdf file containing Ubot & TmBot
%                  that was created from swan2netcdf.m
%   
%   writenetcdf  = accepts either true or false 
%                  if 'True' write tauw and fw to netcdf
%   
%    
%   
%
% OUTPUT:
%    tauw = array containing values of tauw at all time for each node in 
%    model domain
%
% EXAMPLE USAGE
%    tauw = calc_tauw(z0,'skg4.3.nc',true); 
%    
%
% Author(s):  
%    Eric Holmes (University of Massachusetts Dartmouth)
%
% Revision history
%    Initially created 09-24-2010
%   
%==============================================================================

%------------------------------------------------------------------------------
% Convert z0 to meters
%------------------------------------------------------------------------------
z0 = z0/1000;

%------------------------------------------------------------------------------
% Set Ubot & TmBot from netCDF file
%------------------------------------------------------------------------------
nc = netcdf(ncfile,'w');

Ubot = nc{'Ubot'}(:,:);
TmBot = nc{'TmBot'}(:,:);



%------------------------------------------------------------------------------
% Set Generic Values to Salinity and Temperature
% This is temporary, it would be better to use actual salinity and
% temperature data.
%------------------------------------------------------------------------------
vectorsize=size(Ubot);
T=zeros(vectorsize);
S=zeros(vectorsize);

T(:,:)=15;
S(:,:)=30;

%------------------------------------------------------------------------------
% Call Kinematic Viscosity Routine and Density Routine
%------------------------------------------------------------------------------
nu = SW_Kviscosity(T,S); % 
rho = SW_Density(T,S); %


%------------------------------------------------------------------------------
% Calculate fwr & fws according to Soulsby pg. 77-79
%------------------------------------------------------------------------------

%----CONSTANTS-----------------------------------------------------------------
%ks = z0*30.;                    % Nikuradse roughness
A = Ubot.*TmBot/(2.*pi);        % semi-orbital excursion
% r = A/ks;                     % relative roughness
Rw = Ubot.*A./nu;                 % wave Reynolds
%------------------------------------------------------------------------------
fwr = 1.39*(A/z0).^(-0.52);     % case in which flow is rough turbulent 
                                % Soulsby (62a)

                                % Smooth Cases
for i = 1:vectorsize(1,1)
    for j = 1:vectorsize(1,2)
        if(Rw(i,j) > 5e5)
            B=0.0521; N=0.187;  % case in which flow is smooth turbulent
        else
            B=2; N=0.5;         % case in which flow is laminar
        end;
    end
end
                                
fws = B*Rw.^(-N);               % smooth bed friction factor Soulsby (63)

%------------------------------------------------------------------------------
% Choose wave friction factor for current conditions
%------------------------------------------------------------------------------
fw = zeros(vectorsize);

for i = 1:vectorsize(1,1)
    for j = 1:vectorsize(1,2)
        fw(i,j) = max(fwr(i,j),fws(i,j));   
                                % wave friction factor
    end;
end;

tauw = 0.5*rho.*fw.*Ubot.*Ubot;    % wave shear stress
%tauw(isnan(tauw)) = 0;

%------------------------------------------------------------------------------
% Write Values of tauw to ncfile
%------------------------------------------------------------------------------
if (writenetcdf == true)
        
    nc = netcdf(ncfile,'w');
    
    nc{'tauw'} = ncfloat('time','node');
    nc{'tauw'}.long_name = 'Bed Shear Stress';
    nc{'tauw'}.units     = 'N/m^2';
    
    for i=1:vectorsize(1,1)
        nc{'tauw'}(i,1:vectorsize(1,2)) = tauw(i,1:vectorsize(1,2));
    end
    
    close(nc);
    
end;









end
