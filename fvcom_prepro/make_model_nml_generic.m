function make_model_nml_generic(inputConf,yy,mm)
% script to create the model run namelist
% fname='nador_year2046_automatic.nml';
% fnml=fopen(fname,'wt');
% parameters to change for this run
days2sec=24*60*60;

% inputConf.modelYear= 2006;
% inputConf.casename= 'aqua_v14';
% inputConf.base= '/data/medusa/rito/models/FVCOM/runVigo/';
% inputConf.version= 'ver3.2.1'
% inputConf.SIGMA_LEVELS_FILE     = 'sigma_geom.dat';
% inputConf.report= 60 % interval to report to the screen in model minutes
% for mm=8:10;
inputConf.startDate= [yy mm 01 00 00 00 ];
inputConf.endDate=[yy mm+1 01 00 00 00 ];
%% Entries to change
% these cannot be left to the default values
inputConf.EXTSTEP_SECONDS =inputConf.timestep;
inputConf.IRAMP = floor(inputConf.ramp*days2sec./(inputConf.EXTSTEP_SECONDS*inputConf.isplit)); % ramp over one day
inputConf.START_DATE=datestr(inputConf.startDate,'yyyy-mm-dd HH:MM:SS');%           '2046-02-01 00:00:00';
inputConf.END_DATE=datestr(inputConf.endDate,'yyyy-mm-dd HH:MM:SS');% '2046-03-01 00:00:00';
inputConf.RST_FIRST_OUT=inputConf.START_DATE;
% Change sigma file
inputConf.IREPORT = floor(inputConf.report*60./(inputConf.timestep*inputConf.isplit));
inputConf.BOTTOM_ROUGHNESS_FILE  = [inputConf.casename  '_z0=',num2str(inputConf.bed_roughness), '.nc'];
inputConf.INPUT_DIR= inputConf.fvcom_input_nml;
inputConf.OUTPUT_DIR = inputConf.fvcom_output_nml;

inputConf.NC_FIRST_OUT=inputConf.START_DATE;
inputConf.NCAV_FIRST_OUT=inputConf.START_DATE;
inputConf.PROJECTION_REFERENCE=inputConf.projection;
% if isfield(inputConf,'TS_nudge')
%     inputConf.OBC_TEMP_NUDGING_TIMESCALE = 1/(inputConf.TS_nudge*3600/inputConf.timestep);
%     inputConf.OBC_SALT_NUDGING_TIMESCALE = 1/(inputConf.TS_nudge*3600/inputConf.timestep);
%     inputConf.OBC_FABM_NUDGING_TIMESCALE = 1/(inputConf.TS_nudge*3600/inputConf.timestep);
% end



[fmt,nml]=make_default_nml(inputConf);

if isfield(inputConf,lower('OFFLINE_FABM_FILE'))
else
    nml.NML_FABM = rmfield(nml.NML_FABM,'OFFLINE_FABM_FILE');
    
end

% get name of nml blocks
nml_blocks=fieldnames(nml);
nml_vars={};
change_vars=fieldnames(inputConf);
for nn=1:length(nml_blocks)
    nml_vars=fieldnames(nml.(nml_blocks{nn}));
    for vv=1:length(nml_vars)
        var_idx=strcmpi(nml_vars(vv),change_vars);
        if any(var_idx)
            
            change_field=upper(change_vars{find(var_idx)});
            disp(['Changing variable ',nml_blocks{nn},'.',change_field])
            nml.(nml_blocks{nn}).(change_field)=inputConf.(change_vars{find(var_idx)});
        end
    end
end
% % get name of nml blocks
% nml_blocks=fieldnames(nml);
% nml_vars={};
% change_vars=fieldnames(change);
% for nn=1:length(nml_blocks)
%     nml_vars=fieldnames(nml.(nml_blocks{nn}));
%     for vv=1:length(nml_vars)
%         var_idx=strcmp(nml_vars(vv),change_vars);
%         if any(var_idx)
%             disp(['Changing variable ',nml_blocks{nn},'.',change_vars{find(var_idx)}])
%             nml.(nml_blocks{nn}).(change_vars{find(var_idx)})=change.(change_vars{find(var_idx)});
%         end
%     end
% end
res = write_model_nml(inputConf,nml,fmt);
assert(res == 0, 'Error writting namelist file %s', ...
    fullfile(inputConf.fvcom_model, sprintf('%s.nml', inputConf.casename)))
% end
% fprintf(
%  &NML_CASE
%  CASE_TITLE      = 'Nadoor Production FVCOM3.1.6 code'
%  TIMEZONE        = 'UTC';%
%  DATE_FORMAT     = 'YMD'
%  DATE_REFERENCE  = 'default'
%  START_DATE      = '2046-02-01 00:00:00'
%  END_DATE        = '2046-03-01 00:00:00'
% %/
