function [varargout]=read_netCDF_FVCOM(varargin)
%
% Function to extract data from a Netcdf file output from FVCOM
%
% [data] = read_netCDF_FVCOM(varargin)
%
% DESCRIPTION:
%    Function to extract data from a Netcdf file output from FVCOM
%    Outputs data in cell array
%
% INPUT [keyword pairs]:
%   Options are passed in pairs.
%   The list of options (in no particular order) includes:
%   params_opts={'time','data_dir','file_netcdf','varnames','nele_idx','node_idx','siglay_idx','siglev_idx'};
%
%   time  time_interval = {'30/01/06 00:00:00','01/02/06 23:00:00'} or -1 to
%   extract all times in NC file
%
%   data_dir '/home_nfs/rito/models/FVCOM/...' directory where NC file is.
%   default value is ../fvcom_postproc/netcdf
%
%   file_netcdf 'filename.nc'  default value is file_netcdf='*.nc', but it
%   only access the first file in alphabetical order in the directory
%
%   varnames  {'Itime','Itime2','xc','yc','art1','art2','h','siglay','siglev','nv','zeta','ua','va'}
%   cell array of variable names to be read from NC file
%   the variables need to exist in the file but they are case insensitive. The complete list is given by
%   running this script with varnames set to [];
%
%   The variables can be restricted in five possible dimensions
%   node_idx, nele_idx, siglev_idx, siglay_idx and time_idx
%   default values cause the script to extract all available data for all
%   possible dimensions. time_idx is constructed from time_interval. All
%   other indices need to be zero referenced as is Netcdf standard.
%   No checks are done on the bounds of each dimension so make sure you
%   choose them right!
%
%
% OUTPUT:
%    data = cell array with variables in the order they were requested
%
% EXAMPLE USAGE
%   var_2_xtractFVCOM_0 = {'Itime','Itime2','xc','yc','art1','art2','h','siglay','siglev','nv','zeta','ua','va'};
%   date_range={'30/01/06 00:00:00','15/02/06 23:00:00'};  1 if all available data wanted
%   node_idx=[10:30,40:50];% zero referenced!
%   FVCOM_data=read_netCDF_FVCOM_in_progress('time',date_range,'data_dir','/home_nfs/models/FVCOM/runCO2_leak/output/',...
%     'file_netcdf','co2_S1_0001.nc','siglev_idx',1,...
%     'node_idx',node_idx,'varnames',var_2_xtractFVCOM_0);
%
%
% Author(s):
%    Ricardo Torres - Plymouth Marine Laboratory 2012
%    Hakeem Johnson - CH2M
%
% Revision history
%   v0 March 2012
%==============================================================================
%%
%------------------------------------------------------------------------------
%  Parse input arguments
%------------------------------------------------------------------------------
CD=pwd;
disp('Using date conversion of +678942 to go from FVCOM time to matlab time')
time_offset = 678942;
params_opts={'time','data_dir','file_netcdf','varnames','nele_idx','node_idx','siglay_idx','siglev_idx'};

disp('Parameters being used are ...')
var_in_list = {'all_data','netfile_dir','file_netcdf','varnames','nele_idx','node_idx','siglay_idx','siglev_idx'};
all_data = 1;
netfile_dir = '../fvcom_postproc/netcdf';
file_netcdf='*.nc';
siglay_idx=-1;
siglev_idx=-1;
nele_idx=-1;node_idx=-1;
time_idx=-1;
varnames={};
for aa=1:2:nargin
    res=strcmp(varargin(aa),params_opts);
    if ~isempty(res),
        eval([var_in_list{res},' = varargin{aa+1};'])
        disp([params_opts{res}])
    end
end
%------------------------------------------------------------------------------
% sort and remove repeats all indices elements, nodes or layers to increasing values
%------------------------------------------------------------------------------
nele_idx=unique(nele_idx);
node_idx=unique(node_idx);
siglay_idx=unique(siglay_idx);
siglev_idx=unique(siglev_idx);
%
RestrictDims.Name={'node' 'nele' 'siglay' 'siglev' 'time'};
RestrictDims.idx={node_idx, nele_idx, siglay_idx, siglev_idx, time_idx};
%
if ~isempty(varnames)
    nvarnames = length(varnames);
    for nn=1:nvarnames
        data{nn} = [];
    end
end
%%
%------------------------------------------------------------------------------
% Open netcdf file
%------------------------------------------------------------------------------
file_netcdf=[netfile_dir file_netcdf];
filesINdir=dir(file_netcdf);
file_netcdf= fullfile(netfile_dir,filesINdir(1).name);
nc = netcdf.open(file_netcdf, 'NC_NOWRITE');
disp(['NetCDF file ', file_netcdf,' opened successfully.'])
% Get information from netcdf file
info=ncinfo(file_netcdf);
% Extract all possible dimensions in file
DimsAll=info.Dimensions;
% Extract variable names in  nc file
Vars=struct2cell(info.Variables);
vars = squeeze(Vars(1,:,:));
%%
%------------------------------------------------------------------------------
% find variable Itime
%------------------------------------------------------------------------------
try
    Itime.idx=find(strcmpi(vars,'Itime'));
    Itime.ID=netcdf.inqVarID(nc,'Itime');
    Itime.Data  = netcdf.getVar(nc,Itime.ID,'int32');
    Itime2.Data  = netcdf.getVar(nc,Itime.ID+1,'int32');
    %
    [start_d(1),end_d(1)] = deal(double(Itime.Data(1))+time_offset,double(Itime.Data(end))+time_offset);
    [start_d(2),end_d(2)] = deal(double(Itime2.Data(1)),double(Itime2.Data(end)));
    %
    var_time = double(Itime.Data)+time_offset+double(Itime2.Data)./(24*600*6000);
    start_date=sum(start_d.*[1 1/(24*60*60*1000)]);     %hkj missing 1000 inserted
    end_date = sum(end_d.*[1 1/(24*60*60*1000)]);       %hkj missing 1000 inserted
catch me
    fprintf('No ''Itime'' and/or ''Itime2'' variables, using ''time'' instead.\n(%s)\n', me.message)
    Itime.idx=find(strcmpi(vars,'time'));
    Itime.ID=netcdf.inqVarID(nc,'time');
    Itime.Data  = netcdf.getVar(nc,Itime.ID,'double');

    var_time = (Itime.Data)+time_offset;
    [start_date,end_date] = deal(var_time(1),var_time(end));
end

if (length(all_data)==2)
    req_st = datenum(all_data{1},'dd/mm/yy HH:MM:SS');
    req_end = datenum(all_data{2},'dd/mm/yy HH:MM:SS');
else
    req_st =start_date;
    req_end =end_date;
end
time_idx = find(req_st <= var_time &   var_time <= req_end );
% add correct time_idx to RestrictDims
RestrictDims.idx{end}=time_idx;
disp(['Start and end of file, ', datestr(start_date),' ',datestr(end_date)])
%%
%------------------------------------------------------------------------------
% Return information about file to the screen
%------------------------------------------------------------------------------

disp(['Possible variables to extract are: '])
for ii = 1 : length(vars)
    fprintf('%s\n ',vars{ii})
end
if isempty(varnames)
    disp(['Stopping, Choose a variable from the list above : '])
    varargout{1} = 0;
    netcdf.close(nc)
    return
end
%%
%------------------------------------------------------------------------------
% re-organise RestrictDims to follow order of dimensions in nc file from FVCOM
%------------------------------------------------------------------------------
cc=1;
for dd=1:length(DimsAll)
    idx=find(strcmpi(RestrictDims.Name,DimsAll(dd).Name));
    if ~isempty(idx)
        TEMP{cc}=RestrictDims.Name{idx};
        TEMPidx{cc}=RestrictDims.idx{idx};
        cc=cc+1;
    end
end
RestrictDims.Name=TEMP;
RestrictDims.idx=TEMPidx;clear TEMP TEMPidx
%%
%------------------------------------------------------------------------------
% Start Processing extraction of data from NC file
%------------------------------------------------------------------------------
disp(['NetCDF file ', file_netcdf,' opened successfully.'])
nvars=length(info.Variables);
for aa=1:length(varnames)
%------------------------------------------------------------------------------
% Extract number of dimensions, lengths and names of all variables
%------------------------------------------------------------------------------
    TF = strcmpi(varnames{aa},vars);
    if ~isempty(find(TF));
    varidx(aa) = find(TF);TF = sum(TF);
    dimens=ndims(aa);
        disp(['Variable ',vars{varidx(aa)},' found in file'])
    else
        disp(['Variable ',varnames{aa},' NOT found in file Stopping. Check variable names.'])
        netcdf.close(nc)
        varargout{1} = 0;
        return
    end
    varID=netcdf.inqVarID(nc,vars{varidx(aa)});

    [name,xtype,dimids,natts] = netcdf.inqVar(nc,varID);
    dimens=length(dimids);

    for dd=1:length(dimids)
        [dimName{dd},dimLength(dd)] = netcdf.inqDim(nc,dimids(dd));
        disp(['Variable ',name,' has ',num2str(dimens),' dimensions: ',dimName{dd}])
    end
%------------------------------------------------------------------------------
% Get the data!
%------------------------------------------------------------------------------

    start=zeros(size(dimLength));
    count=dimLength;
    switch dimens
        case 1
            % only one dimension present in variable
            switch dimName{1}
                case 'time'
                    if time_idx>=0
                        % only restrict data on access if dimension is TIME
                        %hkj it appears the first value in matlab netcdf interface is 0.
                        %hkj time_idx(1) CORRECTED TO time_idx(1)-1
                        eval([varnames{aa},'=netcdf.getVar(nc,varID,time_idx(1)-1,length(time_idx));'])
                    end
                case 'nele'
                    eval([varnames{aa},'=netcdf.getVar(nc,varID);'])
                    if nele_idx>=0
                        eval([varnames{aa},' = ',varnames{aa},'(nele_idx);'])
                    end
                case 'node'
                    eval([varnames{aa},'=netcdf.getVar(nc,varID);'])
                    if node_idx>=0
                        eval([varnames{aa},' = ',varnames{aa},'(node_idx);'])
                    end
                otherwise
                    disp('Unkown dimension for variable ',name,' Skipping to next one function call');
            end
        otherwise
            % identified dimensions to restrict
            do_restrict=zeros(size(dimName));
            dimidx=nan*ones(size(dimName));
            for dd=1:length(dimName)
                test=find(strcmpi(RestrictDims.Name,dimName{dd}));
                if ~isempty(test); dimidx(dd)=test;   end
            end
            % create start index for dimensions of the variable to
            % access
            if (sum(isfinite(dimidx))==length(dimidx))
                % we have two valid dimension indices, proceed
                for dd=1:length(dimidx)
                    % if restriction is not -1 then select specified
                    % indices otherwise read all
                    if RestrictDims.idx{dimidx(dd)}(1)>=0
                        if (strcmpi(dimName(dd),'time'))
                            start(dd)=RestrictDims.idx{dimidx(dd)}(1)-1;
                            count(dd)=length(start(dd)+1:RestrictDims.idx{dimidx(dd)}(end));
                        else
                            start(dd)=RestrictDims.idx{dimidx(dd)}(1);
                            count(dd)=length(start(dd):RestrictDims.idx{dimidx(dd)}(end));
                        end
                        do_restrict(dd)=1;
                    end
                end
            else
                fprintf('Wrong selection of dimensions to extract, \n  Extracting all values in current variable \n');
            end
            eval([varnames{aa},'=netcdf.getVar(nc,varID,start,count);'])
            % only restrict if required...
            if sum(do_restrict)
                for dd=1:length(do_restrict)
                    sd=dd-1;
                    % calculate indices to extract (might not have been
                    % consecutive numbers)
                    idx=RestrictDims.idx{dimidx(dd)}-start(dd)+1;
                    if ( do_restrict(dd) & ~(count(dd)==length(idx)) )
                        [~,idx]=setdiff(start(dd):RestrictDims.idx{dimidx(dd)}(end),RestrictDims.idx{dimidx(dd)});
                        eval([varnames{aa},' = shiftdim(',varnames{aa},',sd);'])
                        switch  dimens
                            case 2
                                eval([varnames{aa},'(idx, :) = [];'])
                            case 3
                                eval([varnames{aa},'(idx, :,:) = [];'])

                            case 4
                                eval([varnames{aa},'(idx, :,:,:) = [];'])
                        end
                        eval([varnames{aa},' = shiftdim(',varnames{aa},',dimens-sd);'])
                    end
                end
            end
    end
    eval(['data(aa) = {[data{aa};',varnames{aa},']};'])
    eval(['clear ',varnames{aa}])
end
%%
%------------------------------------------------------------------------------
% Tidy up, finish and return data
%------------------------------------------------------------------------------

netcdf.close(nc)
cd(CD)
varargout{1} = data;
return
