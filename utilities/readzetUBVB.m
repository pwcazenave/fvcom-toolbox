% program to read POLCOMS output bin files dailymean
% create file used next in many matlab file
% output matrix  only
% periode 1 for 'janjun' 2 for 'juldec'
% function [u,v,tmp,sal]=readUVT()
% opts = opts.outputfolder;opts.inputfolder;opt.RUN etc
function [data]=readzetUBVB(opts,daysf)
% chose if we had outputed the Vertical velocity, the TKE and the WA of
% polcoms
% variables=1 if they are outputed
% variables=0 otherwise
% outputfolder='/users/modellers/rito/research/Projects/ECOOP/WEC_V0.2/yr2/'
% inputfolder='/users/modellers/rito/research/Projects/ECOOP/WEC_V0.2/yr2/'
readf=opts.zetUVfile;
LEV=opts.PolcomsLevs;
zet=[];u=[];v=[];
%fid is for the input file
readf
fid=fopen(readf,'r','n');
dump= fread(fid,1,'int32');
l= fread(fid,1,'int32')
m= fread(fid,1,'int32')
n= fread(fid,1,'int32')
npsea= fread(fid,1,'int32');
dump = fread(fid,2,'int32');
isea= fread(fid,npsea,'int32');
dump = fread(fid,2,'int32');
jsea= fread(fid,npsea,'int32');
dump = fread(fid,2,'int32');

l= fread(fid,1,'int32');
m= fread(fid,1,'int32');
n= fread(fid,1,'int32');
npusea= fread(fid,1,'int32');
dump = fread(fid,2,'int32');
iusea= fread(fid,npusea,'int32');
dump = fread(fid,2,'int32');
jusea= fread(fid,npusea,'int32');
dump = fread(fid,2,'int32');
hfb = ftell(fid)
for jk=1:daysf
    %         read from  the binary file
    %          hfb_old = ftell(fid)
    u=nan*ones(l,m);
    v=nan*ones(l,m);
    zet=nan*ones(l,m);
    
    itimt4 =fread(fid,1,'int32');disp(['Day :',num2str(jk)]);
    dump = fread(fid,2,'int32');
    dumpzet = fread(fid,[npsea],'float32')';
    dump = fread(fid,2,'int32');
    dumpu = fread(fid,[npusea],'float32')';
    dump = fread(fid,2,'int32');
    dumpv = fread(fid,[npusea],'float32')';
    dump = fread(fid,2,'int32');
    for i=1:npusea
        u(iusea(i),jusea(i))=dumpu(i);
        v(iusea(i),jusea(i))=dumpv(i);
    end
    for i=1:npsea
        zet(isea(i),jsea(i))=dumpzet(i);
    end
    %         print in the matlab file
    if isfield(opts,'PolcomsPoints') && ~isempty(opts.PolcomsPoints)
        
        data.zet(jk,:) = zet(sub2ind(size(v),opts.PolcomsPoints(:,1),opts.PolcomsPoints(:,2)));
        data.ub(jk,:) = u(sub2ind(size(v),opts.PolcomsPoints(:,1),opts.PolcomsPoints(:,2)));
        data.vb(jk,:) = v(sub2ind(size(v),opts.PolcomsPoints(:,1),opts.PolcomsPoints(:,2)));
        
    else
        data.zet(jk,:,:) = zet;
        data.ub(jk,:,:) = u;
        data.vb(jk,:,:) = v;
    end
    %      disp([ftell(fid),ftell(fid)-hfb_old])
end
fclose(fid);
return
