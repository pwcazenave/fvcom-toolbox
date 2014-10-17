function [data,selection] =read_netCDF_FVCOM(varargin)
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
%   possible dimensions. No checks are done on the bounds of each dimension
%   so make sure you choose them right!
%
% OUTPUT:
%    data = struct with fields whose names match those from the list of
%    input variables extracted ('varnames').
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
%   2014-08-20 - Complete the functionality to be able to slice the data
%   along any dimension (siglay, time, node etc.).
%   2014-10-17 - Fix ability to slice with any combination of space
%   (horizontal and vertical) and time.
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
    'node_idx', 'siglay_idx', 'siglev_idx', 'timestride'};

if ftbverbose
    fprintf('Input parameters being used are:\n')
end
var_in_list = {'all_data', 'netfile_dir', 'file_netcdf', 'varnames', ...
    'nele_idx', 'node_idx', 'siglay_idx', 'siglev_idx', 'timestrd'};
all_data = 1;
netfile_dir = '../fvcom_postproc/netcdf';
file_netcdf='*.nc';
siglay_idx=-1;
siglev_idx=-1;
nele_idx=-1;
node_idx=-1;
time_idx=-1;
varnames={};
timestrd=1;
for aa=1:2:nargin
    res=strcmp(varargin(aa),params_opts);
    if sum(res)
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
        data.(varnames{nn}) = [];
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
idx=find(strcmpi(cat(1,{DimsAll.Name}),'time'));
last_entry=DimsAll(idx).Length;
Itime=[];Itime2=[];
% tic
try
    Itime.idx=find(strcmpi(vars,'Itime'));
    Itime.ID=netcdf.inqVarID(nc,'Itime');
    Itime.Data(1)  = netcdf.getVar(nc,Itime.ID,0,1,'int32');
    Itime.Data(2)  = netcdf.getVar(nc,Itime.ID,last_entry-1,1,'int32');
    Itime2.Data(1)  = netcdf.getVar(nc,Itime.ID+1,0,1,'int32');
    Itime2.Data(2)  = netcdf.getVar(nc,Itime.ID+1,last_entry-1,1,'int32');
    
    [start_d(1),end_d(1)] = deal(double(Itime.Data(1))+time_offset,double(Itime.Data(end))+time_offset);
    [start_d(2),end_d(2)] = deal(double(Itime2.Data(1)),double(Itime2.Data(end)));
    
    start_date=sum(start_d.*[1 1/(24*60*60*1000)]);     %hkj missing 1000 inserted
    end_date = sum(end_d.*[1 1/(24*60*60*1000)]);       %hkj missing 1000 inserted
    var_time =  netcdf.getVar(nc,Itime.ID,[0],[10],'double')+time_offset+...
        netcdf.getVar(nc,Itime.ID+1,0,10,'double')./(24*600*6000) ;
    
    DeltaT=median(diff(var_time));
    var_time = start_date:DeltaT:(end_date-DeltaT);
    
catch me
    if ftbverbose
        warning('No ''Itime'' and/or ''Itime2'' variables, using less precise ''time'' instead.\n(%s)\n', me.message)
    end
    Itime.idx=find(strcmpi(vars,'time'));
    Itime.ID=netcdf.inqVarID(nc,'time');
    Itime.Data(1)  = netcdf.getVar(nc,Itime.ID,0,1,'double');
    Itime.Data(2)  = netcdf.getVar(nc,Itime.ID,last_entry-1,1,'double');
    [start_date,end_date] = deal(Itime.Data(1)+time_offset,Itime.Data(end)+time_offset);
    DeltaT=(end_date-start_date)./last_entry;
    var_time = start_date:DeltaT:(end_date-DeltaT);
end
% toc
if length(all_data) == 2
    req_st = datenum(all_data{1},'dd/mm/yy HH:MM:SS');
    req_end = datenum(all_data{2},'dd/mm/yy HH:MM:SS');
else
    req_st = start_date;
    req_end =end_date;
end
time_idx = find(req_st <= var_time &   var_time <= req_end );
time_idx = time_idx(1:timestrd:end);
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
selection=[];
for aa=1:length(varnames)
    selection.(varnames{aa}).start=-1;
    selection.(varnames{aa}).count=-1;
    %----------------------------------------------------------------------
    % Extract number of dimensions, lengths and names of all variables
    %----------------------------------------------------------------------
% tic
    fprintf('Processing variable %s: ', varnames{aa})
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
        error('\nNOT found in file. Stopping. Check input variable names.')
    end
    varID=netcdf.inqVarID(nc,vars{varidx(aa)});

    [name,xtype,dimids,natts] = netcdf.inqVar(nc,varID);
    dimens=length(dimids);

    for dd=1:length(dimids)
        [dimName{dd}, dimLength(dd)] = netcdf.inqDim(nc,dimids(dd));
        if ftbverbose
            if dd == 1
                if length(dimids) == 1
                    fprintf('%i dimension: %s ', dimens, dimName{dd})
                else
                    fprintf('%i dimensions: %s ', dimens, dimName{dd})
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
                        eval([varnames{aa},'=netcdf.getVar(nc,varID,time_idx(1)-1,length(time_idx),timestrd,''double'');'])
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
            dimidx=nan(size(dimName));
            clear start count stride
            for dd=1:length(dimName)
                start.(dimName{dd})=[];
                count.(dimName{dd})=[];
                stride.(dimName{dd})=[];
                test=find(strcmpi(RestrictDims.Name,dimName{dd}));
                if ~isempty(test); dimidx(dd)=test; end
            end
            % create start index for dimensions of the variable to
            % access
            if any(isfinite(dimidx))
                % we have at least two valid dimension indices, proceed
                for dd=1:length(dimidx)
                    % restrict time as range but node and nele dims are
                    % considered as stations rather than ranges.
                    % if restriction is not -1 then select specified
                    % indices otherwise read all
                    if ~isnan(dimidx(dd)) && RestrictDims.idx{dimidx(dd)}(1)>=0
                        if (strcmpi(dimName(dd),'time'))
                            start.(dimName{dd})=RestrictDims.idx{dimidx(dd)}(1)-1;
                            count.(dimName{dd})=length(RestrictDims.idx{dimidx(dd)});
                            stride.(dimName{dd})=timestrd;
                            
                        else
                            for ss=1:length(RestrictDims.idx{dimidx(dd)})
                                start.(dimName{dd})(ss)=RestrictDims.idx{dimidx(dd)}(ss)-1;
                                count.(dimName{dd})(ss)=1;
                                stride.(dimName{dd})=1;
                            end
                        end
                        do_restrict(dd)=1;
                    else
                        start.(dimName{dd})=0;
                        count.(dimName{dd})=dimLength(dd);
                        stride.(dimName{dd})=1;
                    end
                end
            else
                if ftbverbose
                    fprintf('Wrong selection of dimensions to extract.\nExtracting all values in current variable.\n');
                end
            end
            %
            %             eval([varnames{aa},'=netcdf.getVar(nc,varID,start,count,''double'');'])
            cc_names=fieldnames(count);
            clear read_start read_count read_stride
            switch sum(do_restrict) % there are dimensions to restrict
                case 1 % only one dimension to restrict
                    switch find(do_restrict) % find position of restrictive variable
                        case 1 % restrict the first variable
                            % but the variable can have more than 2 dimensions
                            switch dimens
                                % initialise variable
                                case 2
                                    rr=[min(sum(count.(cc_names{1})),dimLength(1)) min(sum(count.(cc_names{2})),dimLength(2))];
                                case 3
                                    rr=[min(sum(count.(cc_names{1})),dimLength(1)),...
                                        min(sum(count.(cc_names{2})),dimLength(2)),...
                                        min(sum(count.(cc_names{3})),dimLength(3))];
                            end

                            eval([varnames{aa},'=nan(rr);'])
                            % reorganize start and count arrays
                            read_start(find(~do_restrict))=start.(cc_names{find(~do_restrict)});
                            read_count(find(~do_restrict))=count.(cc_names{find(~do_restrict)});
                            read_stride(find(~do_restrict))=stride.(cc_names{find(~do_restrict)});
                            
                            for cc=1:length(start.(cc_names{find(do_restrict)}))
                                read_start(find(do_restrict))=start.(cc_names{find(do_restrict)})(cc);
                                read_count(find(do_restrict))=count.(cc_names{find(do_restrict)})(cc);
                                read_stride(find(do_restrict))=stride.(cc_names{find(do_restrict)});
                                
                                var_dump=netcdf.getVar(nc,varID,read_start,read_count,read_stride,'double');
                                
                                eval([varnames{aa},'(cc,:)=var_dump;'])
                                clear var_dump
                            end
                        case 2 % restrict the second variable (ie depth)
                            % but the variable can have more than 2 dimensions
                            switch dimens
                                % initialise variable
                                case 2
                                    rr=[min(sum(count.(cc_names{1})),dimLength(1)) min(sum(count.(cc_names{2})),dimLength(2))];
                                case 3
                                    rr=[min(sum(count.(cc_names{1})),dimLength(1)),...
                                        min(sum(count.(cc_names{2})),dimLength(2)),...
                                        min(sum(count.(cc_names{3})),dimLength(3))];
                            end

                            eval([varnames{aa},'=nan(rr);'])
                            % reorganize start and count arrays
                            read_start(find(~do_restrict))=start.(cc_names{find(~do_restrict)});
                            read_count(find(~do_restrict))=count.(cc_names{find(~do_restrict)});
                            read_stride(find(~do_restrict))=stride.(cc_names{find(~do_restrict)});

                            for cc=1:length(start.(cc_names{logical(do_restrict)}))
                                read_start(find(do_restrict))=start.(cc_names{find(do_restrict)})(cc);
                                read_count(find(do_restrict))=count.(cc_names{find(do_restrict)})(cc);
                                read_stride(find(do_restrict))=stride.(cc_names{find(do_restrict)});
                                var_dump=netcdf.getVar(nc,varID,read_start,read_count,read_stride,'double');
                                try
                                    eval([varnames{aa},'(:,cc)=var_dump;'])
                                catch
                                    eval([varnames{aa},'(:,:,cc)=var_dump;'])

                                end
                                clear var_dump
                            end
                        case 3 % restrict the second variable (ie depth)
                            % but the variable needs to have at least 3 dimensions
                            rr=[min(sum(count.(cc_names{1})),dimLength(1)),...
                                min(sum(count.(cc_names{2})),dimLength(2)),...
                                min(sum(count.(cc_names{3})),dimLength(3))];

                            eval([varnames{aa},'=nan(rr);'])
                            % reorganize start and count arrays
                            % There are now 2 unrestricted dimensions
                            for tt=find(~do_restrict)
                                read_start(tt)=start.(cc_names{tt});
                                read_count(tt)=count.(cc_names{tt});
                                read_stride(tt)=stride.(cc_names{tt});
                                
                            end

                            % check if time is one of them
                            if ~isempty(find(dimidx==5))
                                do_time = find(dimidx==5); % 5 is the index for time
                                % reorganize start and count arrays
                                read_start(do_time)=start.(cc_names{do_time});
                                read_count(do_time)=count.(cc_names{do_time});
                                read_stride(do_time)=stride.(cc_names{do_time});
                                eval([varnames{aa},'=netcdf.getVar(nc,varID,read_start,read_count,read_stride,''double'');'])
                            else % we are looking at stations or depth layers
                                for cc=1:length(start.(cc_names{(do_restrict)}))
                                    read_start(find(do_restrict))=start.(cc_names{find(do_restrict)})(cc);
                                    read_count(find(do_restrict))=count.(cc_names{find(do_restrict)})(cc);
                                    read_stride(find(do_restrict))=stride.(cc_names{find(do_restrict)});
                                    var_dump=netcdf.getVar(nc,varID,read_start,read_count,read_stride,'double');
                                    
                                    switch dimName(find(do_restrict))
                                        case 'node' | 'nele'
                                            eval([varnames{aa},'(cc,:,:)=var_dump;'])
                                        case 'siglay' | 'siglev'
                                            eval([varnames{aa},'(:,cc,:)=var_dump;'])
                                    end
                                    clear var_dump
                                end

                            end
                    end
                    eval(['selection.',varnames{aa},'.start=start;'])
                    eval(['selection.',varnames{aa},'.count=count;'])

                case 2 % Two dimension to restrict!
                    % but the variable can have more than 2 dimensions
                    switch dimens
                        % initialise variable
                        case 2
                            rr=[min(sum(count.(cc_names{1})),dimLength(1)) min(sum(count.(cc_names{2})),dimLength(2))];
                        case 3
                            rr=[min(sum(count.(cc_names{1})),dimLength(1)),...
                                min(sum(count.(cc_names{2})),dimLength(2)),...
                                min(sum(count.(cc_names{3})),dimLength(3))];
                    end

                    eval([varnames{aa},'=nan(rr);'])
                    % check if time is one of them
                    if ~isempty(find(dimidx==5))
                        do_time = find(dimidx==5); % 5 is the index for time
                        % reorganize start and count arrays
                        read_start(do_time)=start.(cc_names{do_time});
                        read_count(do_time)=count.(cc_names{do_time});
                        read_stride(do_time)=stride.(cc_names{do_time});
                        % search for the non_restrictive variable
                        %                         cc=1
                        %                         while ~(length( start.(cc_names{cc}))==1);cc=cc+1;end
                        % esto esta mal.... tengo que incluir otra opcion por si tenemos una
                        % variable de dos dimensiones donde los dos son restrictivas....
                        cc=find(~do_restrict);
                        if isempty(cc);cc=length(cc_names);end
                        read_start(cc)=start.(cc_names{cc});
                        read_count(cc)=count.(cc_names{cc});
                        read_stride(cc)=stride.(cc_names{cc});
                        do_other = setdiff(dimidx,[dimidx(cc),5]) ; % one of these is also restrictive...
                        do_other=find(dimidx==do_other);

                        for cc=1:length(start.(cc_names{do_other}))
                            read_start(do_other)=start.(cc_names{do_other})(cc);
                            read_count(do_other)=count.(cc_names{do_other})(cc);
                            read_stride(do_other)=stride.(cc_names{do_other});
                            var_dump=netcdf.getVar(nc,varID,read_start,read_count,read_stride,'double');
                            switch do_other
                                case 1
                                    eval([varnames{aa},'(cc,:,:)=var_dump;'])
                                case 2
                                    eval([varnames{aa},'(:,cc,:)=var_dump;'])
                                case 3
                                    eval([varnames{aa},'(:,:,cc)=var_dump;'])
                            end
                            clear var_dump
                        end
                    else % time is not one of them so we need to restrict both variables...
                        % in this case it doesn't really matter
                        % which one we restrict firts...

                        for kk=1:length(start.(cc_names{1}))
                            % reorganize start and count arrays
                            read_start(1)=start.(cc_names{1})(kk);
                            read_count(1)=count.(cc_names{1})(kk);
                            read_stride(1)=stride.(cc_names{1});
                            for cc=1:length(start.(cc_names{2}))
                                read_start(2)=start.(cc_names{2})(cc);
                                read_count(2)=count.(cc_names{2})(cc);
                                read_stride(2)=stride.(cc_names{2});
                                var_dump=netcdf.getVar(nc,varID,read_start,read_count,read_stride,'double');
                                
                                eval([varnames{aa},'(kk,cc)=var_dump;'])
                                clear var_dump
                            end

                        end

                    end

                    eval(['selection.',varnames{aa},'.start=start;'])
                    eval(['selection.',varnames{aa},'.count=count;'])
                    
                case 3 % three dimension to restrict!
                    % but the variable can have more than 2 dimensions
                    switch dimens
                        % initialise variable
                        case 2
                            rr=[min(sum(count.(cc_names{1})),dimLength(1)) min(sum(count.(cc_names{2})),dimLength(2))];
                        case 3
                            rr=[min(sum(count.(cc_names{1})),dimLength(1)),...
                                min(sum(count.(cc_names{2})),dimLength(2)),...
                                min(sum(count.(cc_names{3})),dimLength(3))];
                    end
                    
                    eval([varnames{aa},'=nan(rr);'])
                    % check if time is one of them
                    if isempty(find(dimidx==5));disp('This won''t work, try again');return;end
                    do_time = find(dimidx==5); % 5 is the index for time
                    % reorganize start and count arrays
                    read_start(do_time)=start.(cc_names{do_time});
                    read_count(do_time)=count.(cc_names{do_time});
                    read_stride(do_time)=stride.(cc_names{do_time});
                    
                    % search for the non_restrictive variable
                    %                         cc=1
                    %                         while ~(length( start.(cc_names{cc}))==1);cc=cc+1;end
                    % esto esta mal.... tengo que incluir otra opcion por si tenemos una
                    % variable de dos dimensiones donde los dos son restrictivas....
                    [~,do_other] = setdiff(dimidx,[dimidx(do_time)]) ; % these are also restrictive and are not time...
                    if length(count.(cc_names{do_other(1)})) <length(count.(cc_names{do_other(2)}))
                        do_one=do_other(1);do_two=do_other(2);
                    else
                        do_one=do_other(2);do_two=do_other(1);
                    end
                    
                    for cc=1:length(start.(cc_names{do_one}))
                        read_start(do_one)=start.(cc_names{do_one})(cc);
                        read_count(do_one)=count.(cc_names{do_one})(cc);
                        read_stride(do_one)=stride.(cc_names{do_one});
                        
                        for pp=1:length(start.(cc_names{do_two}))
                            read_start(do_two)=start.(cc_names{do_two})(pp);
                            read_count(do_two)=count.(cc_names{do_two})(pp);
                            read_stride(do_two)=stride.(cc_names{do_two});
                            
                            var_dump=netcdf.getVar(nc,varID,read_start,read_count,read_stride,'double');
                            eval([varnames{aa},'(pp,cc,:)=var_dump;'])
                        end
                        clear var_dump
                    end
                    eval(['selection.',varnames{aa},'.start=start;'])
                    eval(['selection.',varnames{aa},'.count=count;'])
                case 0 % there are NO dimensions to restrict and 3 dimensions haven't been coded yet!!
                    
                    for nn=1:length(cc_names)
                        read_start(nn)=start.(cc_names{nn});
                        read_count(nn)=count.(cc_names{nn});
                         read_stride(nn)=stride.(cc_names{nn});
                    end
                    eval([varnames{aa},'=netcdf.getVar(nc,varID,read_start,read_count,read_stride,''double'');'])
                    eval(['selection.',varnames{aa},'.start=start;'])
                    eval(['selection.',varnames{aa},'.count=count;'])

            end
    end
    eval(['data.(varnames{aa}) = ',varnames{aa},';'])
    eval(['clear ',varnames{aa}])
%     toc
end

%--------------------------------------------------------------------------
% Tidy up, finish and return data
%--------------------------------------------------------------------------

netcdf.close(nc)

if ftbverbose
    fprintf('end   : %s \n', subname)
end
