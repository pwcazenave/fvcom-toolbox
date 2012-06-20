function [varargout]=read_netCDF_FVCOM(varargin)
% funtion to read underway NetCDF files from discovery
% possible intputs are
% all_data = 0 single file or 1 for multiple files
% directory where the data is stored
% files = cell array of filenames or if single string then a search pattern
% i.e. *position-4000.gps
% varnames in cell array of variables to extract i.e. lat, lon, time
% they need to be exact matches of variables contained in the netcdf file
% The options could be (but depend on your model output:
% Itime, Itime2, Times, h, iint, kh, km, kq, l, lat, latc, lon, lonc,
% nprocs, nv, partition, q2, q2l, salinity, siglay, siglay_DT,
% siglay_shift_DT, siglev, siglev_DT, temp, time, u, ua, v, va, wet_cells,
% wet_cells_prev_ext, wet_cells_prev_int, wet_nodes, wet_nodes_prev_int, x,
% xc, y, yc, zeta,
% possible attributes for each variable are:
% nc{'measureTS'} = ncdouble('time'); %% 2983 elements.
% nc{'measureTS'}.element_name = ncchar(''measure time'');
% nc{'measureTS'}.cardinalitymin = nclong(1);
% nc{'measureTS'}.cardinalitymax = nclong(1);
% nc{'measureTS'}.comment = ncchar(''time of measure as determined by the GPS'');
% nc{'measureTS'}.long_name = ncchar(''measure timestamp'');
% nc{'measureTS'}.units = ncchar(''day since 1899-12-30T00:00:00 UTC'');
% nc{'measureTS'}.shortunits = ncchar(''days'');
% nc{'measureTS'}.positive = ncchar(''up'');
% nc{'measureTS'}.C_format = ncchar(''%14.7f'');
% nc{'measureTS'}.axis = ncchar(''T'');
% nc{'measureTS'}.measuretimedata = ncchar(''measureTS'');
% nc{'measureTS'}.valid_max = ncdouble(100000);
% nc{'measureTS'}.valid_min = ncdouble(30000);
% nc{'measureTS'}.precision = nclong(12);
% nc{'measureTS'}.scale = nclong(7);
% nc{'measureTS'}.FillValue_ = ncdouble(0);
% nc{'measureTS'}.missing_value = ncdouble(0);
% nc{'measureTS'}.scale_factor = ncdouble(1);
% nc{'measureTS'}.add_offset = ncdouble(0);
% nc{'measureTS'}.element_version = ncchar(''1.0'');
% nc{'measureTS'}.valid_range = ncchar(''30000.000000,100000.000000'');
netcdfpath
addpath D:\Research\Models\UMASSDvisit\matlab;
%%
CD=pwd;
disp(['Using date conversion of +678942 to go from FVCOM time to matlab time'])
time_offset = 678942;
params_opts={'time','data_dir','files','nzopt','trnsopt','varnames','trns_idx','nz_idx'}

disp(['Default values are ...'])
var_in_list = {'all_data','netfile_dir','files','nz','trnsxy','varnames','trns_idx','nz_idx'};
all_data = 1;
netfile_dir = 'D:\research\Data\ICON\underway\netcdf';
files = '*.nc';
nz=0;nz_idx=-1;
trnsxy=0;trns_idx=-1;
varnames={'time','lat','lon','h',};
for aa=1:2:nargin
    res=strmatch(varargin(aa),params_opts)
    if ~isempty(res),
        eval([var_in_list{res},' = varargin{aa+1};'])
        disp([params_opts{res}])
    end
    %     eval([var_in_list{aa},' = varargin{aa};'])
end
if nz==0
    nz_idx=-1;
end
if trnsxy==0;
    trns_idx=-1;
end

nvarnames = length(varnames);
for nn=1:nvarnames
    data{nn} = [];
end
%%
cd (netfile_dir)
filen = dir(files);
for ff=1:length(filen)
    file_netcdf{ff} = filen(ff).name;
end
%     filen(1).name = files;
%     file_netcdf{1} = files;
% Extract time range for all available files
files_to_read = [];times=[];
for ff=1:length(filen)
    nc = netcdf(file_netcdf{ff}, 'nowrite');

    gatts = att(nc); % this are the general file attributes
    vars = var(nc);
    for ii = 1 : length(vars)
        output = name(vars{ii});
        to_show{ii}=output;
        switch output
            case 'Itime'
                [start_d(1),end_d(1)] = deal(vars{ii}(1)+time_offset,vars{ii}(end)+time_offset)
                [start_d(2),end_d(2)] = deal(vars{ii+1}(1),vars{ii+1}(end))
                var_time{ff} = (vars{ii}(:)+time_offset)+(vars{ii+1}(:)./(24*600*6000));
        end
    end
    start_date=sum(start_d.*[1 1/(24*60*60)]);
    end_date = sum(end_d.*[1 1/(24*60*60)]);
    disp(['Start and end of file, ', datestr(start_date),' ',datestr(end_date)])
    files_to_read(ff)=ff;
    times(ff,:) = [start_date ,end_date];
    close(nc)
end
disp(['Possible variables to extract are: '])
for ii = 1 : length(to_show)
    fprintf('%s, ',to_show{ii})
end


if (length(all_data)==2)
    req_st = datenum(all_data{1},'dd/mm/yy HH:MM:SS');
    req_end = datenum(all_data{2},'dd/mm/yy HH:MM:SS');
    files_to_read(:)=NaN;
    sel1 = find(req_st > times(:,1) & req_st < times(:,2));
    sel2 = find(req_end > times(:,1) & req_end < times(:,2));
    sel=[sel1,sel2];
    sel(find(~diff(sel)))=[];
    files_to_read(sel)=1;
    if (req_end > times(end,2));files_to_read(end)=1;end
end
file_netcdf(isnan(files_to_read))=[];
var_time(isnan(files_to_read))=[];
time_idx = find(req_st <= var_time{1} &   var_time{1} <= req_end );

clear vars
for ff=1:length(file_netcdf)
    cdfid = ncmex('OPEN',file_netcdf{ff},'NOWRITE')
    if cdfid ==-1
        disp(['NetCDF file ', file_netcdf{ff},' not found'])
        return
    else
        disp(['NetCDF file ', file_netcdf{ff},' opened successfully.'])
    end
    [ndims, nvars, natts, recdim, status] = ncmex('INQUIRE', cdfid);
    nombre={};dim={};
    for aa=1:nvars
        [nombre{aa}, datatype, ndims(aa), dim{aa}, natts, status] = ncmex('VARINQ', cdfid, aa-1);
    end
    for aa=1:length(varnames)
        TF = strcmpi(varnames{aa},nombre);varidx(aa) = find(TF)-1;TF = sum(TF);
        dimens=ndims(aa);
        if TF;
            disp(['Variable ',varnames{aa},' found in file'])
        else
            disp(['Variable ',varnames{aa},' NOT found in file Stopping. Check variable names.'])
            ncmex('CLOSE',cdfid)
            return
        end
        [dud, dud, dimens, vardims, dud, status] = ncmex('varinq', cdfid, varidx(aa));


        switch dimens
            case 1
                x_dimid = vardims(dimens)

                [dud, x_length, status] = ncmex('diminq', cdfid, x_dimid);
                eval(['[',varnames{aa},', status] = ncmex(''VARGET'', cdfid, varidx(aa), [0], [-1],''autoscale'');'])
            case 2
                z_dimid = vardims(dimens-1)
                x_dimid = vardims(dimens)

                [dud, z_length, status] = ncmex('diminq', cdfid, z_dimid);
                [dud, x_length, status] = ncmex('diminq', cdfid, x_dimid);
                eval(['[',varnames{aa},', status] = ncmex(''VARGET'', cdfid, varidx(aa), [0 0], [-1 -1],''autoscale'');'])
            case 3
                z_dimid = vardims(dimens-1)
                x_dimid = vardims(dimens)
                t_dimid = vardims(dimens-2)

                [dud, z_length, status] = ncmex('diminq', cdfid, z_dimid);
                [dud, x_length, status] = ncmex('diminq', cdfid, x_dimid);
                [dud, t_length, status] = ncmex('diminq', cdfid, t_dimid);
                eval(['[',varnames{aa},', status] = ncmex(''VARGET'', cdfid, varidx(aa), [time_idx(1)-1 nz 0], [length(time_idx) length(nz_idx)*sign(nz_idx(1)) -1],''autoscale'');'])
                if trnsxy

                    eval([varnames{aa},' = squeeze(',varnames{aa},'(trns_idx, : ,:));'])
                end
        end
        eval(['data(aa) = {[data{aa};',varnames{aa},']};'])
        eval(['clear ',varnames{aa}])

    end
    ncmex('CLOSE',cdfid)
end
cd(CD)
varargout{1} = data;
    ncmex('CLOSE',cdfid)
return
%
% pcolor(zero_to_nan(squeeze(value(:,:,1,1)))'),colorbar
%
%
% ncmex('CLOSE',cdfid)
%
% [value, status] = ncmex('VARGET', cdfid, varid, [0 0 0], [londimid latdimid depdimid])
%  [value, status] = ncmex('VARGET', cdfid, varid, [0 0 0], [-1 -1 -1])
%
%
%   status = ncmex('VARPUT1', cdfid, varid, coords, value, autoscale)
%   [value, status] = ncmex('VARGET1', cdfid, varid, coords, autoscale)
%   status = ncmex('VARPUT', cdfid, varid, start, count, value, autoscale)
%   [value, status] = ncmex('VARGET', cdfid, varid, start, count, autoscale)
%   status = ncmex('VARPUTG', cdfid, varid, start, count, stride, [], value, autoscale)
%   [value, status] = ncmex('VARGETG', cdfid, varid, start, count, stride, [], autoscale)
%   status = ncmex('VARRENAME', cdfid, varid, 'name')
%
%   status = ncmex('ATTPUT', cdfid, varid, 'name', datatype, len, value)
%   [datatype, len, status] = ncmex('ATTINQ', cdfid, varid, 'name')
%   [value, status] = ncmex('ATTGET', cdfid, varid, 'name')
%   status = ncmex('ATTCOPY', incdf, invar, 'name', outcdf, outvar)
%   ['name', status] = ncmex('ATTNAME', cdfid, varid, attnum)
%   status = ncmex('ATTRENAME', cdfid, varid, 'name', 'newname')
%   status = ncmex('ATTDEL', cdfid, varid, 'name')
%
%   status = ncmex('RECPUT', cdfid, recnum, [data], autoscale, recdim)
%   [[data], status] = ncmex('RECGET', cdfid, recnum, autoscale, recdim)
%   [[recvarids], [recsizes], status] = ncmex('RECINQ', cdfid, recdim)
%
%   len = ncmex('TYPELEN', datatype)
%   old_fillmode = ncmex('SETFILL', cdfid, fillmode)
%
%   old_ncopts = ncmex('SETOPTS', ncopts)
%   ncerr = ncmex('ERR')
%   code = ncmex('PARAMETER', 'NC_...')
%
%   Notes:
%    1. The rcode is always zero.
%    2. The dimid can be number or name.
%    3. The varid can be number or name.
%    4. The attname can be name or number.
%    5. The operation and parameter names are not case-sensitive.
%    6. The cmode defaults to 'NC_NOCLOBBER'.
%    7. The mode defaults to 'NC_NOWRITE'.
%    8. The value -1 determines length automatically.
%    9. The operation names can prepend 'nc'.
%   10. The parameter names can drop 'NC_' prefix.
%   11. Dimensions: Matlab (i, j, ...) <==> [..., j, i] NetCDF.
%   12. Indices and identifiers are zero-based.
%   13. One-dimensional arrays are returned as column-vectors.
