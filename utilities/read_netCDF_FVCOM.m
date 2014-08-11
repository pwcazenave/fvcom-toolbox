function data =read_netCDF_FVCOM(varargin)
% Function to extract data from a Netcdf file output from FVCOM.
%
% data = read_netCDF_FVCOM(varargin)
%
% DESCRIPTION:
%    Function to extract data from a netCDF file output from FVCOM. Outputs
%    data in cell array.
%
% INPUT [keyword pairs]:
%   Options are passed in pairs.
%
%   The list of keywords is:
%       - 'time'
%       - 'data_dir'
%       - 'file_netcdf'
%       - 'varnames'
%       - 'nele_idx'
%       - 'node_idx'
%       - 'siglay_idx'
%       - 'siglev_idx'
%
%   'time' - {'30/01/06 00:00:00', '01/02/06 23:00:00'} or -1 to extract
%   all times.
%
%   'data_dir' - '/home/fvcom/results/...' directory where netCDF file is.
%   Default value is ../fvcom_postproc/netcdf
%
%   'file_netcdf' - 'filename.nc'. Default value is '*.nc', but it only
%   access the first file in alphabetical order in the directory.
%
%   'varnames' - Cell array of variable names to read from the netCDF file.
%   The variables need to exist in the file but they are case insensitive.
%   Choose FVCOM output variable names. For example:
%       - 'Itime'
%       - 'Itime2'
%       - 'xc'
%       - 'yc'
%       - 'art1'
%       - 'art2'
%       - 'h'
%       - 'siglay'
%       - 'siglev'
%       - 'nv'
%       - 'zeta'
%       - 'ua'
%       - 'va'
%   The complete list for a given file is given by running this script with
%   varnames set to [].
%
%   The variables can be restricted in five possible dimensions:
%       - 'node_idx'
%       - 'nele_idx'
%       - 'siglev_idx'
%       - 'siglay_idx'
%       - 'time_idx'
%   Default values cause the script to extract all available data for all
%   possible dimensions. All indices except time need to be zero referenced
%   as is standard for netCDF indexing. No checks are done on the bounds of
%   each dimension so make sure you choose them right!
%
%
% OUTPUT:
%    data = cell with variables in the order they were requested.
%
% EXAMPLE USAGE
%   vars = {'Times', 'xc', 'yc', 'h', 'siglay', 'nv', 'zeta', 'ua', 'va'};
%   date_range = {'30/01/06 00:00:00', '15/02/06 23:00:00'};
%   node_idx = [10:30, 40:50]; % zero referenced!
%   data_dir = '/home/fvcom/results/output/';
%   FVCOM = read_netCDF_FVCOM('data_dir', data_dir, ...
%       'file_netcdf', 'casename_0001.nc', ...
%       'time', date_range, ...
%       'siglev_idx', 1, ...
%       'node_idx', node_idx, ...
%       'varnames', vars);
%
% Author(s):
%    Ricardo Torres - Plymouth Marine Laboratory 2012
%    Hakeem Johnson - CH2M
%    Pierre Cazenave - Plymouth Marine Laboratory
%
% Revision history:
%   v0 March 2012
%   2014-03-06 - Add the global verbose flag. Also tidy up the help a bit.
%   Also change some verbose statements to use fprintf instead of disp for
%   better control over formatting. Also fixed a bug where if a 2D array
%   was requested after a 3D array, the 2D array would cause the function
%   to crash (because it was using a 3D index for getVar).
%
%==========================================================================

global ftbverbose
subname = 'read_netCDF_FVCOM';

if ftbverbose
    fprintf('\nbegin : %s \n', subname)
end

%--------------------------------------------------------------------------
%  Parse input arguments
%--------------------------------------------------------------------------

params_opts = {'time', 'data_dir', 'file_netcdf', 'varnames', 'nele_idx', ...
    'node_idx', 'siglay_idx', 'siglev_idx'};

if ftbverbose
    fprintf('Input parameters being used are:\n')
end
var_in_list = {'all_data', 'netfile_dir', 'file_netcdf', 'varnames', ...
    'nele_idx', 'node_idx', 'siglay_idx', 'siglev_idx'};
all_data = 1;
netfile_dir = '../fvcom_postproc/netcdf';
file_netcdf='*.nc';
siglay_idx=-1;
siglev_idx=-1;
nele_idx=-1;
node_idx=-1;
time_idx=-1;
varnames={};
for aa=1:2:nargin
    res=strcmp(varargin(aa),params_opts);
    if ~isempty(res)
        eval([var_in_list{res},' = varargin{aa+1};'])
        if ftbverbose
            fprintf(' %s\n', params_opts{res})
        end
    end
end

%--------------------------------------------------------------------------
% Sort (and remove repeats) for all indices elements, nodes or layers
%--------------------------------------------------------------------------
nele_idx=unique(nele_idx);
node_idx=unique(node_idx);
siglay_idx=unique(siglay_idx);
siglev_idx=unique(siglev_idx);

RestrictDims.Name={'node' 'nele' 'siglay' 'siglev' 'time'};
RestrictDims.idx={node_idx, nele_idx, siglay_idx, siglev_idx, time_idx};

if ~isempty(varnames)
    nvarnames = length(varnames);
    for nn=1:nvarnames
        data{nn} = [];
    end
end

%--------------------------------------------------------------------------
% Open netcdf file
%--------------------------------------------------------------------------
file_netcdf=fullfile(netfile_dir, file_netcdf);
filesINdir=dir(file_netcdf);
file_netcdf= fullfile(netfile_dir,filesINdir(1).name);
nc = netcdf.open(file_netcdf, 'NC_NOWRITE');
if ftbverbose
    if length(file_netcdf) > 50
        % Truncate output file name to display.
        fprintf('NetCDF file ...%s opened successfully.\n', file_netcdf(end-70:end))
    else
        fprintf('NetCDF file %s opened successfully.\n', file_netcdf)
    end
end
% Get information from netcdf file
info=ncinfo(file_netcdf);
% Extract all possible dimensions in file
DimsAll=info.Dimensions;
% Extract variable names in  nc file
Vars=struct2cell(info.Variables);
vars = squeeze(Vars(1,:,:));

%--------------------------------------------------------------------------
% Find variable Itime
%--------------------------------------------------------------------------
if ftbverbose
    fprintf('Using date conversion of +678942 days to go from FVCOM time (Modified Julian Day) to MATLAB time.\n')
end
time_offset = 678942;

try
    Itime.idx=find(strcmpi(vars,'Itime'));
    Itime.ID=netcdf.inqVarID(nc,'Itime');
    Itime.Data  = netcdf.getVar(nc,Itime.ID,'int32');
    Itime2.Data  = netcdf.getVar(nc,Itime.ID+1,'int32');

    [start_d(1),end_d(1)] = deal(double(Itime.Data(1))+time_offset,double(Itime.Data(end))+time_offset);
    [start_d(2),end_d(2)] = deal(double(Itime2.Data(1)),double(Itime2.Data(end)));

    var_time = double(Itime.Data)+time_offset+double(Itime2.Data)./(24*600*6000);
    start_date=sum(start_d.*[1 1/(24*60*60*1000)]);     %hkj missing 1000 inserted
    end_date = sum(end_d.*[1 1/(24*60*60*1000)]);       %hkj missing 1000 inserted
catch me
    if ftbverbose
        warning('No ''Itime'' and/or ''Itime2'' variables, using less precise ''time'' instead.\n(%s)\n', me.message)
    end
    Itime.idx=find(strcmpi(vars,'time'));
    Itime.ID=netcdf.inqVarID(nc,'time');
    Itime.Data  = netcdf.getVar(nc,Itime.ID,'double');

    var_time = (Itime.Data)+time_offset;
    [start_date,end_date] = deal(var_time(1),var_time(end));
end

if length(all_data) == 2
    req_st = datenum(all_data{1},'dd/mm/yy HH:MM:SS');
    req_end = datenum(all_data{2},'dd/mm/yy HH:MM:SS');
else
    req_st = start_date;
    req_end =end_date;
end
time_idx = find(req_st <= var_time &   var_time <= req_end );
% Add correct time_idx to RestrictDims
RestrictDims.idx{end}=time_idx;
if ftbverbose
    fprintf('Start and end of file: %s - %s\n', datestr(start_date), datestr(end_date))
end

%--------------------------------------------------------------------------
% Return information about file to the screen
%--------------------------------------------------------------------------

if ftbverbose
    fprintf('Possible variables to extract are:\n')
end
for ii = 1:length(vars)
    if ftbverbose
        fprintf(' %s\n', vars{ii})
    end
end
if isempty(varnames)
    data = 0;
    netcdf.close(nc)
    error('Stopping. Choose a variable from the list above.')
end

%--------------------------------------------------------------------------
% Re-organise RestrictDims to follow order of dimensions in nc file from
% FVCOM
%--------------------------------------------------------------------------
cc=1;
for dd=1:length(DimsAll)
    idx=find(strcmpi(RestrictDims.Name,DimsAll(dd).Name));
    if ~isempty(idx)
        TEMP{cc}=RestrictDims.Name{idx};
        TEMPidx{cc}=RestrictDims.idx{idx};
        cc=cc+1;
    end
end
RestrictDims.Name = TEMP;
RestrictDims.idx = TEMPidx;
clear TEMP TEMPidx

%--------------------------------------------------------------------------
% Start Processing extraction of data from NC file
%--------------------------------------------------------------------------

for aa=1:length(varnames)
    %----------------------------------------------------------------------
    % Extract number of dimensions, lengths and names of all variables
    %----------------------------------------------------------------------

    % Tidy up the previous iteration's variables so we don't get confused.
    clear dimName dimLength

    TF = strcmpi(varnames{aa},vars);
    if ~isempty(find(TF));
        varidx(aa) = find(TF);
        TF = sum(TF);
        dimens=ndims(aa);
        if ftbverbose
            fprintf('Variable %s found', vars{varidx(aa)})
        end
    else
        netcdf.close(nc)
        varargout{1} = 0;
        error('Variable %s NOT found in file. Stopping. Check input variable names.\n', varnames{aa})
    end
    varID=netcdf.inqVarID(nc,vars{varidx(aa)});

    [name,xtype,dimids,natts] = netcdf.inqVar(nc,varID);
    dimens=length(dimids);

    for dd=1:length(dimids)
        [dimName{dd}, dimLength(dd)] = netcdf.inqDim(nc,dimids(dd));
        if ftbverbose
            if dd == 1
                if length(dimids) == 1
                    fprintf(' with %i dimension: %s ', dimens, dimName{dd})
                else
                    fprintf(' with %i dimensions: %s ', dimens, dimName{dd})
                end
            else
                fprintf('%s ', dimName{dd})
            end
        end
    end
    if ftbverbose; fprintf('\n'); end

    %----------------------------------------------------------------------
    % Get the data!
    %----------------------------------------------------------------------

    start=zeros(size(dimLength));
    count=dimLength;
    switch dimens
        case 1
            % only one dimension present in variable
            switch dimName{1}
                case 'time'
                    if time_idx>=0
                        % Only restrict data on access if dimension is TIME

                        % hkj it appears the first value in matlab netcdf
                        % interface is 0.
                        % hkj time_idx(1) CORRECTED TO time_idx(1)-1.
                        eval([varnames{aa},'=netcdf.getVar(nc,varID,time_idx(1)-1,length(time_idx),''double'');'])
                    end
                case 'nele'
                    eval([varnames{aa},'=netcdf.getVar(nc,varID,''double'');'])
                    if nele_idx>=0
                        eval([varnames{aa},' = ',varnames{aa},'(nele_idx);'])
                    end
                case 'node'
                    eval([varnames{aa},'=netcdf.getVar(nc,varID,''double'');'])
                    if node_idx>=0
                        eval([varnames{aa},' = ',varnames{aa},'(node_idx);'])
                    end
                otherwise
                    if ftbverbose
                        fprintf('Unkown dimension for variable %s. Skipping to next one in function call.\n', name);
                    end
            end
        otherwise
            % identified dimensions to restrict
            do_restrict=zeros(size(dimName));
            dimidx=nan*ones(size(dimName));
            for dd=1:length(dimName)
                test=find(strcmpi(RestrictDims.Name,dimName{dd}));
                if ~isempty(test); dimidx(dd)=test; end
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
                if ftbverbose
                    fprintf('Wrong selection of dimensions to extract.\nExtracting all values in current variable.\n');
                end
            end

            if strcmpi(varnames{aa}, 'Times')
                % A string variable, so don't convert to double.
                eval([varnames{aa},'=netcdf.getVar(nc,varID,start,count);'])
            else
                eval([varnames{aa},'=netcdf.getVar(nc,varID,start,count,''double'');'])
            end

            % only restrict if required...
            if sum(do_restrict)
                for dd=1:length(do_restrict)
                    sd=dd-1;
                    % calculate indices to extract (might not have been
                    % consecutive numbers)
                    idx=RestrictDims.idx{dimidx(dd)}-start(dd)+1;
                    if (do_restrict(dd) && ~(count(dd)==length(idx)))
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

%--------------------------------------------------------------------------
% Tidy up, finish and return data
%--------------------------------------------------------------------------

netcdf.close(nc)

if ftbverbose
    fprintf('end   : %s \n', subname)
end