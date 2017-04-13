function write_FVCOM_tsobcERSEM(tsOBCFile,FileExists,time,nSiglay,Mobj,ERSEMdata,varargin)
% Export temperature and salinity forcing at the open boundary.
%
% function write_FVCOM_tsobc(basename,time,nSiglay,in_temp,in_salt)
%
% DESCRIPTION:
%    Setup an FVCOM hydrographic open boundary forcing file. Supply either
%    uniform values for temperature and salinity or 3D arrays (node,
%    sigma_layers, time).
%
% INPUT
%    tsOBCFile - Filename to use as output
%    FileExists - True if tsOBCFile already exists and we only want to add
%    ERSEM variables to it.
%    time - Time (Modified Julian Days)
%    nSiglay - Number of sigma layers
%    Mobj - Mesh Object with the following fields:
%       - nObs - number of open boundaries
%       - read_obc_nodes - open boundary node cell array (length = nObs)
%       - siglay - sigma layer definitions
%       - siglev - sigma level definitions
%    ERSEMdata - struct with information to output to Netcdf File

%    Optional keyword-argument pairs:
%
%    'strtime' = set to true to output the 'Times' variable
%    'inttime' = set to true to output the 'Itime' and 'Itime2' variables
%    'floattime' = set to true to output the 'time' variable
%
%    This script defaults to writing 'Times' only.
%
%    FVCOM needs only one of:
%        1. Times: character string of times
%        2. Itime and Itime2: integer days and milliseconds since midnight
%        3. time: float days.
%    FVCOM checks for these in the order above and this script defaults to
%    writing Times only. Adjust the keyword-argument pairs to your liking.
%
% OUTPUT:
%    FVCOM hydrographic open boundary netCDF file
%
% Author(s):
%    Geoff Cowles (University of Massachusetts Dartmouth)
%    Pierre Cazenave (Plymouth Marine Laboratory)
%    Karen Amoudry (National Oceanography Centre, Liverpool)
%    Ricardo Torres (Plymouth Marine Laboratory)
% PWC Revision history
%    2012-06-15 Added support for native MATLAB NetCDF routines. Requires
%    MATLAB 2010a or higher.
%    2012-07-16 Removed hard-coded nSiglay and nSiglev and instead moved to
%    arguments list.
%    2012-10-08 Updated help to reflect the fact nSiglev is calculated as
%    nSiglay+1.
%    2012-11-09 Add new arguments to use user defined temperature and
%    salinity.
%    2013-01-09 Add support for 3D input temperature and salinity (such as
%    might be generated with get_POLCOMS_tsobc.m.
%    2016-05-25 Removed the reads of the ASCII configuration files (which
%    was very slow) and instead extracted the relevant information from the
%    supplied mesh object. As such, the requirements for the mesh object
%    have changed, so hopefully this won't bite too many people in the
%    behind. Also simplified the allocation of the arrays when uniform
%    values are given (i.e. when in_salt and in_temp are scalars).
%
% KJA Revision history
%    Undated - Add better check for the size of the input arrays (works
%    with scalars).
%    2013-08-16 - Updated output of Itime2 to avoid rounding errors
%    when converting from double to single format.
% RJT Revision history
%    2013-12-05 Added functionality to output ERSEM nutrients
%    2017-04-10 Option to write to an existing file... i.e. a Nesting file
%==============================================================================

[~, subname] = fileparts(mfilename('fullpath'));
global ftbverbose
if ftbverbose
    fprintf('\nbegin : %s\n', subname)
end
NNuts=length(ERSEMdata)
% Default to string times as FVCOM looks for these first.
strtime = true;
inttime = false;
floattime = false;
for vv = 1:2:length(varargin)
    switch varargin{vv}
        case 'strtime'
            strtime = true;
        case 'inttime'
            inttime = true;
        case 'floattime'
            floattime = true;
    end
end

obc_nodes = [Mobj.read_obc_nodes{:}]';
obc_h = Mobj.h(obc_nodes);
siglev = Mobj.siglev(obc_nodes, :);
siglay = Mobj.siglay(obc_nodes, :);

nTimes = length(time);
nObc = length(obc_nodes);
nSiglev = nSiglay + 1;

if(ftbverbose); fprintf('obc reading complete\n');end;

%--------------------------------------------------------------------------
% Generate the requisite data
%--------------------------------------------------------------------------


% We need to make sigma level and layer data resolved for each node on the
% open boundary (in case we have hybrid coordinates).
if isvector(siglev)
    siglev = repmat(siglev, [nObc, 1]);
end
if isvector(siglay)
    siglay = repmat(siglay, [nObc, 1]);
end

% Check we've got everything the right size and shape.
if nSiglev ~= size(Mobj.(ERSEMdata(1).name), 2) + 1 || size(siglev, 2) ~= size(Mobj.(ERSEMdata(1).name), 2) + 1 || size(siglev, 2) ~= size(Mobj.(ERSEMdata(2).name), 2) + 1
    error('Specified number sigma levels does not match supplied data')
end
if nSiglay ~= size(Mobj.(ERSEMdata(1).name), 2) || size(siglay, 2) ~= size(Mobj.(ERSEMdata(1).name), 2) || size(siglay, 2) ~= size(Mobj.(ERSEMdata(2).name), 2)
    error('Specified number of sigma layers does not match supplied data')
end


%--------------------------------------------------------------------------
% Set netCDF variables and dump to file
%--------------------------------------------------------------------------
if FileExist
    % open boundary forcing
    nc = netcdf.open(tsOBCFile, 'WRITE');
    % read dimensions from the
    % define dimensions
    dimids = netcdf.inqDimIDs(nc);
    for dd =1:length(dimids)
        dimidname=netcdf.inqDim(nc,dimids(dd));
        switch dimidname
            case{'node'}
                
                nobc_dimid=netcdf.inqDimID(nc,dimidname);
            case{'DateStrLen'}
                datestrlen_dimid=netcdf.inqDimID(nc,dimidname);
            case{'time'}
                
                time_dimid=netcdf.inqDimID(nc,dimidname);
                % read time in case is different from data in file
                file_times = netcdf.getVar(nc,netcdf.inqVarID(nc,'Times'))';
                file_timem = datenum(file_times,'yyyy-mm-dd HH:MM:SS.FFF');
            case{'siglay'}
                
                siglay_dimid=netcdf.inqDimID(nc,dimidname);
            case{'siglev'}
                siglev_dimid=netcdf.inqDimID(nc,dimidname);
        end
    end
    netcdf.reDef(nc)
else
    % open boundary forcing
    nc = netcdf.create(tsOBCFile, 'clobber');
    
    % define global attributes
    netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'title','Open boundary ERSEM nudging')
    netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'type','FVCOM TIME SERIES OBC FABM FILE')
    netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'history', sprintf('File created with %s from the MATLAB fvcom-toolbox', subname))
    
    
    % define dimensions
    nobc_dimid=netcdf.defDim(nc,'nobc',nObc);
    datestrlen_dimid=netcdf.defDim(nc,'DateStrLen',26);
    time_dimid=netcdf.defDim(nc,'time',netcdf.getConstant('NC_UNLIMITED'));
    siglay_dimid=netcdf.defDim(nc,'siglay',nSiglay);
    siglev_dimid=netcdf.defDim(nc,'siglev',nSiglev);
    
    % variables
    if strtime
        Times_varid=netcdf.defVar(nc,'Times','NC_CHAR',[datestrlen_dimid, time_dimid]);
        netcdf.putAtt(nc,Times_varid,'time_zone','UTC');
    end
    
    if floattime
        time_varid=netcdf.defVar(nc,'time','NC_FLOAT',time_dimid);
        netcdf.putAtt(nc,time_varid,'long_name','time');
        netcdf.putAtt(nc,time_varid,'units','days since 1858-11-17 00:00:00');
        netcdf.putAtt(nc,time_varid,'time_zone','UTC');
    end
    
    if inttime
        itime_varid=netcdf.defVar(nc,'Itime','NC_INT',time_dimid);
        netcdf.putAtt(nc,itime_varid,'units','days since 1858-11-17 00:00:00');
        netcdf.putAtt(nc,itime_varid,'format','modified julian day (MJD)');
        netcdf.putAtt(nc,itime_varid,'time_zone','UTC');
        
        itime2_varid=netcdf.defVar(nc,'Itime2','NC_INT',time_dimid);
        netcdf.putAtt(nc,itime2_varid,'units','msec since 00:00:00');
        netcdf.putAtt(nc,itime2_varid,'time_zone','UTC');
    end
    
    nobc_varid=netcdf.defVar(nc,'obc_nodes','NC_INT',nobc_dimid);
    netcdf.putAtt(nc,nobc_varid,'long_name','Open Boundary Node Number');
    netcdf.putAtt(nc,nobc_varid,'grid','obc_grid');
    netcdf.putAtt(nc,nobc_varid,'type','data');
    
    obc_h_varid=netcdf.defVar(nc,'obc_h','NC_FLOAT',nobc_dimid);
    netcdf.putAtt(nc,obc_h_varid,'long_name','Open Boundary Depth');
    netcdf.putAtt(nc,obc_h_varid,'units','m');
    netcdf.putAtt(nc,obc_h_varid,'grid','obc_grid');
    netcdf.putAtt(nc,obc_h_varid,'type','data');
    
    obc_siglev_varid=netcdf.defVar(nc,'siglev','NC_FLOAT',[nobc_dimid,siglev_dimid]);
    netcdf.putAtt(nc,obc_siglev_varid,'long_name','ocean_sigma/general_coordinate');
    netcdf.putAtt(nc,obc_siglev_varid,'grid','obc_grid');
    
    obc_siglay_varid=netcdf.defVar(nc,'siglay','NC_FLOAT',[nobc_dimid,siglay_dimid]);
    netcdf.putAtt(nc,obc_siglay_varid,'long_name','ocean_sigma/general_coordinate');
    netcdf.putAtt(nc,obc_siglay_varid,'grid','obc_grid');
end


% nutrients here
for nuts=1:NNuts
    varidN{nuts}=['obc_',ERSEMdata(nuts).name,'_varid'];
    eval([varidN{nuts},'=netcdf.defVar(nc,''',ERSEMdata(nuts).name,''',''NC_FLOAT'',[nobc_dimid,siglay_dimid,time_dimid]);'])
    eval(['netcdf.putAtt(nc,',varidN{nuts},',''long_name'',''',ERSEMdata(nuts).long_name,''');'])
    eval(['netcdf.putAtt(nc,',varidN{nuts},',''units'',''',ERSEMdata(nuts).units,''');'])
    eval(['netcdf.putAtt(nc,',varidN{nuts},',''grid'',''obc_grid'');'])
    % obc_salinity_varid=netcdf.defVar(nc,'obc_salinity','NC_FLOAT',[nobc_dimid,siglay_dimid,time_dimid]);
    % netcdf.putAtt(nc,obc_salinity_varid,'long_name','sea_water_salinity');
    % netcdf.putAtt(nc,obc_salinity_varid,'units','PSU');
    % netcdf.putAtt(nc,obc_salinity_varid,'grid','obc_grid');
end





% end definitions
netcdf.endDef(nc);

if ~FileExists
    % grid and time information already exist in file...
    % write data
    netcdf.putVar(nc,nobc_varid,obc_nodes);
    netcdf.putVar(nc,obc_h_varid,obc_h);
    netcdf.putVar(nc,obc_siglev_varid,siglev);
    netcdf.putVar(nc,obc_siglay_varid,siglay);
    if strtime
        nStringOut = char();
        [nYr, nMon, nDay, nHour, nMin, nSec] = mjulian2greg(time);
        for i=1:nTimes
            nDate = [nYr(i), nMon(i), nDay(i), nHour(i), nMin(i), nSec(i)];
            nStringOut = [nStringOut, sprintf('%04i/%02i/%02i %02i:%02i:%09.6f', nDate)];
        end
        netcdf.putVar(nc,Times_varid,[0, 0],[26, nTimes],nStringOut);
    end
    if floattime
        netcdf.putVar(nc,time_varid,0,numel(time),time);
    end
    if inttime
        netcdf.putVar(nc,itime_varid,floor(time));
        % netcdf.putVar(nc,itime2_varid,0,numel(time),mod(time,1)*24*3600*1000); % PWC original
        % KJA edit: avoids rounding errors when converting from double to single
        % Rounds to nearest multiple of the number of msecs in an hour
        netcdf.putVar(nc,itime2_varid,0,numel(time),round((mod(time,1)*24*3600*1000)/(3600*1000))*(3600*1000));
    end
    for nuts=1:NNuts
        eval(['netcdf.putVar(nc,',varidN{nuts},',Mobj.(ERSEMdata(nuts).name));'])
    end
    
else % file exist and time could be different... check and interpolate if necessary
    if length(time) < length(file_timem)
        for nuts=1:NNuts
            data= Mobj.(ERSEMdata(nuts).name);
            dataint = nans(size(data, 1),size(data, 2),length(file_timem));
            for nn=1:size(data, 1)
                [X1,Y1]=meshgrid(file_timem,siglay(nn,:));
                [X,Y]=meshgrid(time,siglay(nn,:));
                % interpolate ERSEMdata...
                dataint(nn,:,:) = interp2(X,Y,squeeze(data(nn,:,:),X1,Y1));
            end
            
            netcdf.putVar(nc,varidN{nuts},dataint);
        end
    else
        % everything is in the same time frequency
        for nuts=1:NNuts
            eval(['netcdf.putVar(nc,',varidN{nuts},',Mobj.(ERSEMdata(nuts).name));'])
        end
        
    end
end
% close file
netcdf.close(nc);

if ftbverbose
    fprintf('end   : %s\n', subname)
end
