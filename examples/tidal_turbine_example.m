%
% Example code showing how to make a tidal turbine parameter input netCDF
% file for FVCOM.
%
% Rory O'Hara Murray
% 6 Oct 2015
%
clear, close all

% load in PFOW model grid and bathymetry
Mobj = read_fvcom_mesh(['\\isilonml\Shelf_Model\PFOW\1_HD Model\4_Climatology Results\input\PFOW_SMS_1_grd.dat']);
Mobj.h =  read_fvcom_bath('\\isilonml\Shelf_Model\PFOW\1_HD Model\4_Climatology Results\input\PFOW_SMS_1_dep.dat');

lonc = nodes2elems(Mobj.x,Mobj); % calculate approximate x-coordinate at the centre of each element
latc = nodes2elems(Mobj.y,Mobj); % calculate approximate y-coordinate at the centre of each element
[xc yc] = ll2utm(lonc, latc, 30);

% calculate depth at the centre of each element
Mobj.hc = mean(Mobj.h(Mobj.tri),2);

% setup empty turbine array
turbine.numbers     = single(zeros(Mobj.nElems,1));
turbine.sigma_layer = single(zeros(Mobj.nElems,10));
turbine.area        = single(zeros(Mobj.nElems,1));
% turbine.thrust      = single(zeros(Mobj.nElems,1));
% turbine.drag        = single(zeros(Mobj.nElems,1));

% Example turbine positions in the Inner Sound of Stroma
TurbinePositionsll = [-3.137,58.661; ...
-3.12945337127092,58.6596303572796; ...
-3.12960283663035,58.6600269324843; ...
-3.13193517522843,58.658758560404; ...
-3.1320846649985,58.659155132772];

[TurbinePositions(:,1) TurbinePositions(:,2)] = ll2utm(TurbinePositionsll(:,1), TurbinePositionsll(:,2), 30);

% For each turbine find which element it is in and count the number of
% turbines in each element
I = fun_nearest2D(TurbinePositions(:,1), TurbinePositions(:,2), xc, yc);
for ii=1:length(I)
    turbine.numbers(I(ii)) = turbine.numbers(I(ii)) + 1;
end
turbinesI = turbine.numbers>0;
numTurbineElems = sum(turbinesI);       % total number of elements with turbines in
numTurbines = sum(turbine.numbers);     % total number of turbines


% loop through each turbine and work out the fraction spread over each sigma layer
for ii=find(turbinesI)'
    turbine.sigma_layer(ii,:) = turbine_area_sigma(Mobj.hc(ii), 15, 10);
end

% an alternative (simpler) way of specifying the fractional split across
% sigma layers
% turbine.sigma_layer(turbinesI,:) = ones(numTurbineElems,1)*[0 0 0 0 0 0.25 0.25 0.25 0.25 0];

% turbine rotor swept area - 10 m radius.
turbine.area(turbinesI) = 10*10*pi;

% The thrust curve is currently not inputted into FVCOM, this could be
% added in the future though.
% turbine.thrust_curve = [0.99 1 1.25 1.5 1.75 2 2.25 2.5 2.75 3 3.25 3.5 3.75 4 4.01; 0 0.85 0.85 0.85 0.85 0.85 0.85 0.85 0.635 0.490 0.385 0.308 0.250 0.205 0];

% write the netcdf input file
write_FVCOM_TT(turbine,['Tidal_Turbines_Example.nc'],'Example Scenario with 5 tidal turbines in the inner sound');

%% plot the location of the turbines
%plot_fvcom_field(Mobj,elems2nodes(turbine.numbers,Mobj), 'grd', 'w') 
plot_fvcom_field(Mobj,single(turbinesI), 'grd', 'w') 
colormap( jet(2) );
daspect([1 1 1])


