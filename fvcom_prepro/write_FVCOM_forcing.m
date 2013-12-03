function write_FVCOM_forcing(Mobj, fileprefix, data, infos, fver)
% Write data out to FVCOM NetCDF forcing file.
%
% write_FVCOM_forcing(Mobj, fvcom_forcing_file, data, infos, fver)
%
% DESCRIPTION:
%   Takes the given interpolated data (e.g. from grid2fvcom) and writes out
%   to a NetCDF file.
%
% INPUT:
%   Mobj - MATLAB mesh object containing fields:
%       tri - triangulation table for the unstructured grid
%       nVerts - number of grid vertices (nodes)
%       nElems - number of grid elements
%       nativeCoords - model coordinate type ('cartesian' or 'spherical')
%       x, y or lon, lat - node positions (depending on nativeCoords value)
%   fileprefix - Output NetCDF file prefix (plus path) will be
%       fileprefix_{wnd,hfx,evap}.nc if fver is '3.1.0', otherwise output
%       files will be fileprefix_wnd.nc.
%   data - Struct of the data to be written out.
%   infos - Additional remarks to be written to the "infos" NetCDF variable
%   fver - Output for version 3.1.0 or 3.1.6. The latter means all the
%       forcing can go in a single file, the former needs separate files
%       for specific forcing data (wind, heating and precipitation).
%
% The fields in data may be called any of:
%     - 'u10', 'v10', 'uwnd', 'vwnd' - wind components
%     - 'slp'       - sea level pressure
%     - 'pevpr'     - evaporation
%     - 'prate'     - precipitation
%     - 'slp'       - sea level pressure
%     - 'nlwrs'     - net longwave radiation*,**
%     - 'nswrs'     - net shortwave radiation*,**
%     - 'shtfl'     - sensible heat net flux*,**
%     - 'lhtfl'     - latent heat net flux*,**
%     - 'lon'       - longitude (vector)
%     - 'lat'       - latitude (vector)
%     - 'x'         - eastings (vector)
%     - 'y'         - northings (vector)
%     - 'nshf'      - pre-computed net surface heat flux**
%
% Fields marked with an * are combined to form the "surface net heat flux"
% (nshf) as follows:
%
%   nshf = -nlwrs + -nswrs - lhtfl - shtfl;
%
% ** Alternatively, a new field 'nshf' (net surface heat flux) can be
% supplied, in which case shtfl and lhtfl are not necessary as their only
% function is in the calculation of net surface heat flux. This approach
% eliminates the need to interpolate so many variables onto the FVCOM grid,
% thus decreasing the time needed to generate the FVCOM input files.
%
% OUTPUT:
%   FVCOM wind speed forcing NetCDF file(s)
%
% EXAMPLE USAGE:
%   windBase = '/path/to/output/casename_wnd.nc';
%   write_FVCOM_forcing(Mobj, windBase, data, ...
%       'FVCOM atmospheric forcing data', '3.1.6');
%
% Author(s):
%   Pierre Cazenave (Plymouth Marine Laboratory)
%   Karen Thurston (National Oceanography Centre, Liverpool)
%
% PWC Revision history:
%   2012-11-01 - First version based on the parts of grid2fvcom_U10V10.m
%   which dealt with writing the NetCDF file. This version now dynamically
%   deals with varying numbers of forcing data.
%   2012-11-09 - Add the correct calculation for the net surface heat flux.
%   2013-02-14 - Fix the surface net heat flux sign convention, include the
%   longwave radiation variable and add support for a new field in the
%   input struct ('nshf') which contains the pre-computed net surface heat
%   flux.
%   2013-05-13 - Fix the evaporation to use the correct variable from NCEP
%   (pevpr rather than P_E which is actually the precipitation minus the
%   evaporation in Et). The data in Et are calcaulated from lhtfl whereas
%   pevpr comes directly from NCEP and to me it seems more sensible to use
%   that to maintain consistency.
%   2013-05-14 - Add example usage to the help and specify which fields are
%   required in Mobj.
%
% KJT Revision history:
%   2013-01-16 - Added support for output of sea level pressure.
%
% ROM Revision History:
%   Change the output of tri' to tri, as tri was being written the wrong
%   way around (for windows matlab anyway)
%==========================================================================

multi_out = false; % default to 3.1.6, single output file
if nargin < 4 || nargin > 5
    error('Incorrect number of arguments')
elseif nargin == 5
    if strcmpi(fver, '3.1.0')
        multi_out = true;
    end
end

subname = 'write_FVCOM_forcing';

global ftbverbose;
if ftbverbose
	fprintf('\nbegin : %s \n', subname)
end

tri = Mobj.tri;
nNodes = Mobj.nVerts;
nElems = Mobj.nElems;
ntimes = numel(data.time);

if strcmpi(Mobj.nativeCoords, 'cartesian')
    x = Mobj.x;
    y = Mobj.y;
else
    x = Mobj.lon;
    y = Mobj.lat;
end
% Create a string for each variable's coordinate attribute
coordString = sprintf('FVCOM %s coordinates', Mobj.nativeCoords);

% Create element vertices positions
xc = nodes2elems(x, Mobj);
yc = nodes2elems(y, Mobj);

%--------------------------------------------------------------------------
% Create the NetCDF header for the FVCOM forcing file
%--------------------------------------------------------------------------

if multi_out
    suffixes = {'_wnd', '_hfx', '_evap', '_air_press'};
else
    suffixes = {'_wnd'};
end

% We use this variable to indicate whether we were given precalculated net
% surface heat flux.
nshf = 0;

for i=1:length(suffixes)
    nc = netcdf.create([fileprefix, suffixes{i}, '.nc'], 'clobber');

%     netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'type','FVCOM Forcing File')
    netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'title','FVCOM Forcing File')
    netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'institution','Plymouth Marine Laboratory')
    netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'source','FVCOM grid (unstructured) surface forcing')
    netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'history',['File created on ', datestr(now, 'yyyy-mm-dd HH:MM:SS'), ' with write_FVCOM_forcing.m from the MATLAB fvcom-toolbox (https://github.com/pwcazenave/fvcom-toolbox)'])
    netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'references','http://fvcom.smast.umassd.edu, http://codfish.smast.umassd.edu')
    netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'Conventions','CF-1.0')
%     netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'infos',infos)
    netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'CoordinateSystem',Mobj.nativeCoords)
    netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'CoordinateProjection','init=epsg:4326') % WGS84?

    % Dimensions
    nele_dimid=netcdf.defDim(nc,'nele',nElems);
    node_dimid=netcdf.defDim(nc,'node',nNodes);
    three_dimid=netcdf.defDim(nc,'three',3);
    time_dimid=netcdf.defDim(nc,'time',netcdf.getConstant('NC_UNLIMITED'));
    datestrlen_dimid=netcdf.defDim(nc,'DateStrLen',26);

    % Space variables
    x_varid=netcdf.defVar(nc,'x','NC_FLOAT',node_dimid);
    netcdf.putAtt(nc,x_varid,'long_name','nodal x-coordinate');
    netcdf.putAtt(nc,x_varid,'units','meters');

    y_varid=netcdf.defVar(nc,'y','NC_FLOAT',node_dimid);
    netcdf.putAtt(nc,y_varid,'long_name','nodal y-coordinate');
    netcdf.putAtt(nc,y_varid,'units','meters');

    xc_varid=netcdf.defVar(nc,'xc','NC_FLOAT',nele_dimid);
    netcdf.putAtt(nc,xc_varid,'long_name','zonal x-coordinate');
    netcdf.putAtt(nc,xc_varid,'units','meters');

    yc_varid=netcdf.defVar(nc,'yc','NC_FLOAT',nele_dimid);
    netcdf.putAtt(nc,yc_varid,'long_name','zonal y-coordinate');
    netcdf.putAtt(nc,yc_varid,'units','meters');

    nv_varid=netcdf.defVar(nc,'nv','NC_INT',[nele_dimid, three_dimid]);
    netcdf.putAtt(nc,nv_varid,'long_name','nodes surrounding element');

    % Time variables
    time_varid=netcdf.defVar(nc,'time','NC_FLOAT',time_dimid);
    netcdf.putAtt(nc,time_varid,'long_name','time');
    netcdf.putAtt(nc,time_varid,'units','days since 1858-11-17 00:00:00');
    netcdf.putAtt(nc,time_varid,'format','modified julian day (MJD)');
    netcdf.putAtt(nc,time_varid,'time_zone','UTC');

    itime_varid=netcdf.defVar(nc,'Itime','NC_INT',time_dimid);
    netcdf.putAtt(nc,itime_varid,'units','days since 1858-11-17 00:00:00');
    netcdf.putAtt(nc,itime_varid,'format','modified julian day (MJD)');
    netcdf.putAtt(nc,itime_varid,'time_zone','UTC');

    itime2_varid=netcdf.defVar(nc,'Itime2','NC_INT',time_dimid);
    netcdf.putAtt(nc,itime2_varid,'units','msec since 00:00:00');
    netcdf.putAtt(nc,itime2_varid,'time_zone','UTC');
    netcdf.putAtt(nc,itime2_varid,'long_name','time');

    times_varid=netcdf.defVar(nc,'Times','NC_CHAR',[datestrlen_dimid,time_dimid]);
	netcdf.putAtt(nc,times_varid,'long_name','Calendar Date');
    netcdf.putAtt(nc,times_varid,'format','String: Calendar Time');
    netcdf.putAtt(nc,times_varid,'time_zone','UTC');

    % Since we have a dynamic number of variables in the struct, try to be a
    % bit clever about how to create the output variables.
    fnames = fieldnames(data);
    used_varids = cell(0);
    used_fnames = cell(0);
    used_dims = cell(0); % exclude time (assume all variables vary in time)

    for vv=1:length(fnames)
        % Need to check both whether we have the data but also whether
        % we're outputting to several NetCDF files. If so, we drop some
        % variables if we're in the wrong file loop.
        switch fnames{vv}
            case {'uwnd', 'u10'}
                if strcmpi(suffixes{i}, '_wnd') || ~multi_out
                    % wind components (assume we have v if we have u)

                    % On the elements
                    u10_varid=netcdf.defVar(nc,'U10','NC_FLOAT',[nele_dimid, time_dimid]);
                    netcdf.putAtt(nc,u10_varid,'long_name','Eastward Wind Speed');
%                     netcdf.putAtt(nc,u10_varid,'standard_name','Eastward Wind Speed');
                    netcdf.putAtt(nc,u10_varid,'units','m/s');
                    netcdf.putAtt(nc,u10_varid,'grid','fvcom_grid');
                    netcdf.putAtt(nc,u10_varid,'coordinates',coordString);
                    netcdf.putAtt(nc,u10_varid,'type','data');

                    v10_varid=netcdf.defVar(nc,'V10','NC_FLOAT',[nele_dimid, time_dimid]);
                    netcdf.putAtt(nc,v10_varid,'long_name','Northward Wind Speed');
%                     netcdf.putAtt(nc,v10_varid,'standard_name','Northward Wind Speed');
                    netcdf.putAtt(nc,v10_varid,'units','m/s');
                    netcdf.putAtt(nc,v10_varid,'grid','fvcom_grid');
                    netcdf.putAtt(nc,v10_varid,'coordinates',coordString);
                    netcdf.putAtt(nc,v10_varid,'type','data');

                    uwind_varid=netcdf.defVar(nc,'uwind_speed','NC_FLOAT',[nele_dimid, time_dimid]);
                    netcdf.putAtt(nc,uwind_varid,'long_name','Eastward Wind Speed');
                    netcdf.putAtt(nc,uwind_varid,'standard_name','Wind Speed');
                    netcdf.putAtt(nc,uwind_varid,'units','m/s');
                    netcdf.putAtt(nc,uwind_varid,'grid','fvcom_grid');
                    netcdf.putAtt(nc,uwind_varid,'type','data');

                    vwind_varid=netcdf.defVar(nc,'vwind_speed','NC_FLOAT',[nele_dimid, time_dimid]);
                    netcdf.putAtt(nc,vwind_varid,'long_name','Northward Wind Speed');
                    netcdf.putAtt(nc,vwind_varid,'standard_name','Wind Speed');
                    netcdf.putAtt(nc,vwind_varid,'units','m/s');
                    netcdf.putAtt(nc,vwind_varid,'grid','fvcom_grid');
                    netcdf.putAtt(nc,vwind_varid,'type','data');

                    % Only on the elements (both U10/V10 and uwind_speed and
                    % vwind_speed).
                    used_varids = [used_varids, {'u10_varid', 'v10_varid', 'uwind_varid', 'vwind_varid'}];
                    used_fnames = [used_fnames, {'uwnd', 'vwnd', 'uwnd', 'vwnd'}];
                    used_dims = [used_dims, {'nElems', 'nElems', 'nElems', 'nElems'}];
                end
                
            case {'vwnd', 'v10'}
                % We dealt with these in the u component section above, so
                % just pass silently.
                true;

            case 'slp'
                if strcmpi(suffixes{i}, '_air_press') || ~multi_out
                    % Sea level pressure
                    slp_varid=netcdf.defVar(nc,'air_pressure','NC_FLOAT',[node_dimid, time_dimid]);
                    netcdf.putAtt(nc,slp_varid,'long_name','Surface air pressure');
                    netcdf.putAtt(nc,slp_varid,'units','Pa');
                    netcdf.putAtt(nc,slp_varid,'grid','fvcom_grid');
                    netcdf.putAtt(nc,slp_varid,'coordinates',coordString);
                    netcdf.putAtt(nc,slp_varid,'type','data');

                    used_varids = [used_varids, 'slp_varid'];
                    used_fnames = [used_fnames, fnames{vv}];
                    used_dims = [used_dims, 'nNodes'];
                end

            case {'pevpr', 'Et'}
                if strcmpi(suffixes{i}, '_evap') || ~multi_out
                    % Evaporation
                    pevpr_varid=netcdf.defVar(nc,'evap','NC_FLOAT',[node_dimid, time_dimid]);
                    netcdf.putAtt(nc,pevpr_varid,'long_name','Evaporation');
                    netcdf.putAtt(nc,pevpr_varid,'description','Evaporation, ocean lose water is negative');
                    netcdf.putAtt(nc,pevpr_varid,'units','m s-1');
                    netcdf.putAtt(nc,pevpr_varid,'grid','fvcom_grid');
                    netcdf.putAtt(nc,pevpr_varid,'coordinates',coordString);
                    netcdf.putAtt(nc,pevpr_varid,'type','data');

                    used_varids = [used_varids, 'pevpr_varid'];
                    used_fnames = [used_fnames, fnames{vv}];
                    used_dims = [used_dims, 'nNodes'];
                end

            case {'prate', 'P_E'}
                if strcmpi(suffixes{i}, '_evap') || ~multi_out
                    % Precipitation
                    prate_varid=netcdf.defVar(nc,'precip','NC_FLOAT',[node_dimid, time_dimid]);
                    netcdf.putAtt(nc,prate_varid,'long_name','Precipitation');
                    netcdf.putAtt(nc,prate_varid,'description','Precipitation, ocean lose water is negative');
                    netcdf.putAtt(nc,prate_varid,'units','m s-1');
                    netcdf.putAtt(nc,prate_varid,'grid','fvcom_grid');
                    netcdf.putAtt(nc,prate_varid,'coordinates',coordString);
                    netcdf.putAtt(nc,prate_varid,'type','data');

                    used_varids = [used_varids, 'prate_varid'];
                    used_fnames = [used_fnames, fnames{vv}];
                    used_dims = [used_dims, 'nNodes'];
                end

            case 'nswrs'
                if strcmpi(suffixes{i}, '_hfx') || ~multi_out
                    % Shortwave radiation
                    nswrs_varid=netcdf.defVar(nc,'short_wave','NC_FLOAT',[node_dimid, time_dimid]);
                    netcdf.putAtt(nc,nswrs_varid,'long_name','Short Wave Radiation');
                    netcdf.putAtt(nc,nswrs_varid,'units','W m-2');
                    netcdf.putAtt(nc,nswrs_varid,'grid','fvcom_grid');
                    netcdf.putAtt(nc,nswrs_varid,'coordinates',coordString);
                    netcdf.putAtt(nc,nswrs_varid,'type','data');

                    used_varids = [used_varids, 'nswrs_varid'];
                    used_fnames = [used_fnames, fnames{vv}];
                    used_dims = [used_dims, 'nNodes'];
                end

            case 'nlwrs'
                if strcmpi(suffixes{i}, '_hfx') || ~multi_out
                    % Longwave radiation
                    nlwrs_varid=netcdf.defVar(nc,'long_wave','NC_FLOAT',[node_dimid, time_dimid]);
                    netcdf.putAtt(nc,nlwrs_varid,'long_name','Long Wave Radiation');
                    netcdf.putAtt(nc,nlwrs_varid,'units','W m-2');
                    netcdf.putAtt(nc,nlwrs_varid,'grid','fvcom_grid');
                    netcdf.putAtt(nc,nlwrs_varid,'coordinates',coordString);
                    netcdf.putAtt(nc,nlwrs_varid,'type','data');

                    used_varids = [used_varids, 'nlwrs_varid'];
                    used_fnames = [used_fnames, fnames{vv}];
                    used_dims = [used_dims, 'nNodes'];
                end

            case {'shtfl', 'lhtfl', 'nshf'} % , 'nlwrs', 'nswrs'}
                % We can't trigger on nlwrs and nswrs here because they're
                % the triggers for the net longwave and shortwave variables
                % above. Instead, we need to use the latent heat flux
                % variables as the triggers for the net heat flux.
                % We also need to check for the existence of the "net
                % surface heat flux ('nshf')" field which can be created before
                % calling grid2fvcom. This approach means there's fewer
                % calls to the (expensive) interpolation as the net surface
                % heat flux is calculated before being interpolated onto
                % the FVCOM grid. Set the nshf variable accordingly.
                if strcmpi(fnames{vv}, 'nshf')
                    nshf = 1;
                end
                try
                    % We might have already made this attribute, so fail
                    % elegantly if we do. This is because we need to put
                    % all four of shtfl, lhtfl, nlwrs and nswrs to make
                    % Surface Net Heat Flux.
                    if strcmpi(suffixes{i}, '_hfx') || ~multi_out
                        % Surface net heat flux
                        nhf_varid=netcdf.defVar(nc,'net_heat_flux','NC_FLOAT',[node_dimid, time_dimid]);
                        netcdf.putAtt(nc,nhf_varid,'long_name','Surface Net Heat Flux');
                        netcdf.putAtt(nc,nhf_varid,'units','W m-2');
                        netcdf.putAtt(nc,nhf_varid,'grid','fvcom_grid');
                        netcdf.putAtt(nc,nhf_varid,'coordinates',coordString);
                        netcdf.putAtt(nc,nhf_varid,'type','data');
                    end
                end
                if strcmpi(suffixes{i}, '_hfx') || ~multi_out
                    % We need to save the current variable name even if we've
                    % already made its attribute.
                    used_varids = [used_varids, 'nhf_varid'];
                    used_fnames = [used_fnames, fnames{vv}];
                    used_dims = [used_dims, 'nNodes'];
                end

            case {'time', 'lon', 'lat', 'x', 'y'}
                continue

            otherwise
                if ftbverbose
                    warning('Unknown or possibly unused input data type: %s', fnames{vv})
                end
        end
    end

    % End definitions
    netcdf.endDef(nc);

    % Put the easy ones in first.
    netcdf.putVar(nc,nv_varid, tri);
    netcdf.putVar(nc,time_varid,0,ntimes,data.time);
    netcdf.putVar(nc,itime_varid,0,ntimes,floor(data.time));
    netcdf.putVar(nc,itime2_varid,0,ntimes,mod(data.time,1)*24*3600*1000);
    netcdf.putVar(nc,x_varid,x);
    netcdf.putVar(nc,y_varid,y);
    netcdf.putVar(nc,xc_varid,xc);
    netcdf.putVar(nc,yc_varid,yc);

    % Build the time string and output to NetCDF.
    nStringOut = char();
    for tt=1:ntimes
        [nYr, nMon, nDay, nHour, nMin, nSec] = mjulian2greg(data.time(tt));
        nDate = [nYr, nMon, nDay, nHour, nMin, nSec];
        nStringOut = [nStringOut, sprintf('%04i/%02i/%02i %02i:%02i:%02i       ',nDate)];
    end
    netcdf.putVar(nc,times_varid,nStringOut);

    % Now do the dynamic ones. Set the heat flux to not done (hf_done = 0)
    % until we hit one of the holy quad (shtfl, lhtfl, nlwrs and nswrs).
    hf_done = 0;
    for ff=1:length(used_fnames)
        if ftbverbose
            fprintf('write : %s... ', used_fnames{ff})
        end
        if strcmpi(used_fnames{ff}, 'shtfl') || strcmpi(used_fnames{ff}, 'lhtfl') || strcmpi(used_fnames{ff}, 'nlwrs') || strcmpi(used_fnames{ff}, 'nswrs')

            hf_done = hf_done + 1;

            if hf_done == 4 && nshf == 0
                if ftbverbose
                    fprintf('combining heat flux ... ')
                end
                % We've got all four heat parameters, so dump them into the
                % file. We have to flip the signs of the net fluxes and
                % subtract the latent and sensible heat fluxes because
                % NCEP's convention is positive upwards (out of the ocean)
                % and negative downwards (into the ocean), whereas FVCOM is
                % positive downwards (into the ocean) and negative upwards
                % (out of the ocean).
                hf = -data.nlwrs.node + -data.nswrs.node - data.lhtfl.node - data.shtfl.node;
                %hf = data.shtfl.node + data.lhtfl.node + data.nlwrs.node + data.nswrs.node;
                netcdf.putVar(nc,nhf_varid,[0,0],[nNodes,ntimes],hf)
            elseif strcmpi(used_fnames{ff}, 'nswrs') || strcmpi(used_fnames{ff}, 'nlwrs')
                % We've already done the net surface heat flux but we're on
                % either of the other fluxes (short/long wave) which we
                % need to dump. Do that here.
                if strcmpi(used_dims{ff}, 'nNodes')
                    eval(['netcdf.putVar(nc,',used_varids{ff},',[0,0],[',used_dims{ff},',ntimes],data.',used_fnames{ff},'.node);'])
                else
                    eval(['netcdf.putVar(nc,',used_varids{ff},',[0,0],[',used_dims{ff},',ntimes],data.',used_fnames{ff},'.data);'])
                end
            else
                % We haven't got the precomputed net surface heat flux but
                % we haven't yet got enough of the parameters to export the
                % heat flux from the short + long + latent + sensible.
                % Essentially this loop just does hf_done = hf_done + 1.
            end
        elseif strcmpi(used_fnames{ff}, 'nshf') && nshf == 1
            if ftbverbose
                fprintf('existing combined heat flux ... ')
            end
            % We have pre-computed net surface heat flux, in which case set
            % hf_done to 4 and put the data into the netCDF. Also set the
            % nshf variable 1 to stop the net surface heat flux variable
            % being overwritten above.
            hf_done = 4;
            netcdf.putVar(nc, nhf_varid, [0, 0], [nNodes, ntimes], data.nshf.node)
        else
            % One of the other data sets for which we can simply dump the
            % existing array without waiting for other data
            if strcmpi(used_dims{ff}, 'nNodes')
                eval(['netcdf.putVar(nc,',used_varids{ff},',[0,0],[',used_dims{ff},',ntimes],data.',used_fnames{ff},'.node);'])
            else
                eval(['netcdf.putVar(nc,',used_varids{ff},',[0,0],[',used_dims{ff},',ntimes],data.',used_fnames{ff},'.data);'])
            end
        end
        if ftbverbose
            fprintf('done.\n')
        end
    end
    if hf_done < 4
        % hf_done might be higher than four, but unless it is at least
        % four, we haven't got everything we need.
        warning('Did not have all the required heat flux parameters. Need ''shtfl'', ''lhtfl'', ''nlwrs'' and ''nwsrs''.')
    end

    % Close the NetCDF file(s)
    netcdf.close(nc);
end

if ftbverbose
    fprintf('end   : %s \n', subname)
end
