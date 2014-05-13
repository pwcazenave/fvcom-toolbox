function params=read_polcoms_params(fileparams)
% function to read polcoms parameters
% fileparams='/users/modellers/rito/Models/MEDINA/polcoms/MEDI29.parameters'
% list of variables to keep
varlist={'l' 'm' 'n' 'daldi' 'dbedi' 'along1' 'alat1'}
fid = fopen(fileparams);
C = textscan(fid, '%s%s%*[^\n]');
fclose(fid);
for vv=1:length(varlist)
    iv=find(strcmpi(C{2},varlist(vv)));
    params.(varlist{vv})=str2num(C{1}{iv});
end

return