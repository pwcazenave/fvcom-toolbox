function [polcoms]=calc_scoord(polcoms)
% %n=18;iesub=63;jesub=36;
%dir_path='/users/modellers/rito/research/Projects/MERSEA/dataWCH/WCH2_rito/data2003/';
%filenameb= 'WCH.bathy_WCH2'
%filenameiu= 'WCH.ipexu_WCH2'
%filenameib= 'WCH.ipexb_WCH2'
% program to calculate the scoordinates from model run
%dir_data = '/data/milkyway/rito/lolo/mpi_H_POLCOMS_VIGO_3D_OK_2/setups/data2002_03/'
% Model parameters
%n=18;iesub=50;jesub=70;
%  S coordinate parameters
hs = polcoms.bathy;
% Read ipexb

polcoms.ipexbM= zeros(polcoms.iesub,polcoms.jesub);
for aa=1:length(polcoms.isea)
    polcoms.ipexbM(polcoms.isea(aa),polcoms.jsea(aa))=1;
end
polcoms.ipexuM= zeros(polcoms.iesub,polcoms.jesub);
for aa=1:length(polcoms.npusea)
    polcoms.ipexuM(polcoms.iusea(aa),polcoms.jusea(aa))=1;
end

%  Evenly spaced S values
%

dsc = 1.0d0/(polcoms.params.n-2);
sval(1)   = -1.0d0;
for k=2:polcoms.params.n-2
    sval(k) = sval(k-1)+dsc;
end
sval(polcoms.params.n-1) = 0.0d0;
%JPO  START
sval(polcoms.params.n)   = 0.0d0;
%JPO  END
for j=1:polcoms.jesub
    for i=1:polcoms.iesub
        sigo(1,i,j) = -1.0d0;
        for k=2:polcoms.params.n-2
%             c               for identical to sigma on shelf
%             c
            if (hs(i,j)>polcoms.scoord.hc) 
                ffh=(hs(i,j)-polcoms.scoord.hc)/hs(i,j);
                cs = (1.0d0-polcoms.scoord.bb)*(sinh(polcoms.scoord.theta*sval(k)))/sinh(polcoms.scoord.theta)+...
                    polcoms.scoord.bb*(tanh(polcoms.scoord.theta*(sval(k)+0.5d0))...
                    -tanh(0.5d0*polcoms.scoord.theta))/...
                    (2.*tanh(0.5d0*polcoms.scoord.theta));
                sigo(k,i,j) = sval(k)+ffh*(cs-sval(k));
            else
                sigo(k,i,j) = sval(k);
            end
        end
        sigo(polcoms.params.n-1,i,j) = 0.0d0;
%         CJPO START
        sigo(polcoms.params.n,i,j) = 0.0d0;
%         CJPO END
        %
    end
    %
end
%
%     Define coordinates on B points
%
ds = 0.0;
for j=1:polcoms.jesub
    for i=1:polcoms.iesub
        for k=1:polcoms.params.n-2
            ds(k,i,j) = sigo(k+1,i,j)-sigo(k,i,j);
        end
        %c
        %c           Set surface level to ensure correct sum for ds
        %c
        sds = 0.0d0;
        for k=1:polcoms.params.n-3
            sds = sds+ds(k,i,j);
        end
        ds(polcoms.params.n-2,i,j) = 1.0d0-sds;
        %c
        dsu(1,i,j) = ds(1,i,j);
        for k=2:polcoms.params.n-2
            dsu(k,i,j) = 0.5d0*(ds(k,i,j)+ds(k-1,i,j));
        end
        dsu(polcoms.params.n-1,i,j) = ds(polcoms.params.n-2,i,j);
        %c
        %sig(0,i,j) = -1.0d0;
        sig(1,i,j) = -1.0d0+0.5d0*ds(1,i,j);
        for k=2:polcoms.params.n-2
            sig(k,i,j) = sig(k-1,i,j)+dsu(k,i,j);
        end
        sig(polcoms.params.n-1,i,j) = 0.0d0;
        %CJPO START
        sig(polcoms.params.n,i,j) = 0.0d0;
        %CJPO END
    end
end
sigov = nan*ones(size(sig));
dsv = nan*ones(size(dsu));
dsuv = nan*ones(size(dsu));
%c
%c     Average coordinates onto U points
%c
for j=2:polcoms.jesub
    for i=2:polcoms.iesub
        if ( polcoms.ipexbM(i,j)~=0 | polcoms.ipexuM(i,j)~=0 )
            sigov(1,i,j) = -1.0d0;
            for k=2:polcoms.params.n-2
                sigov(k,i,j) = 0.25d0*...
                    (sigo(k,i  ,j  )+sigo(k,i-1,j  )...
                    +sigo(k,i-1,j-1)+sigo(k,i  ,j-1));
            end
            sigov(polcoms.params.n-1,i,j) = 0.0d0;
            %c
        end
        %c
        for k=1:polcoms.params.n-2
            dsv(k,i,j) = sigov(k+1,i,j)-sigov(k,i,j);
        end
        %c
        %cjth        Set surface level to ensure correct sum for DS
        %c
        sds=0.;
        for k=1:polcoms.params.n-3
            sds=sds+dsv(k,i,j);
        end
        dsv(polcoms.params.n-2,i,j) = 1.0d0-sds;
        %c
        dsuv(1,i,j) = dsv(1,i,j);
        dsuv(polcoms.params.n-1,i,j) = dsv(polcoms.params.n-2,i,j);
        for k=2:polcoms.params.n-2
            dsuv(k,i,j) = 0.5d0*(dsv(k,i,j)+dsv(k-1,i,j));
        end
        %c
%        sigv(0,i,j) = -1.0d0;
        sigv(1,i,j) = -1.0d0+0.5d0*dsv(1,i,j);
        for k=2:polcoms.params.n-2
            sigv(k,i,j) = sigv(k-1,i,j)+dsuv(k,i,j);
        end
        sigv(polcoms.params.n-1,i,j) = 0.0d0;
        %CJPO START
        sigv(polcoms.params.n,i,j) = 0.0d0;
        %CJPO END
    end
end
sigv(:,1,:)=sigo(:,1,:);sigv(:,:,1)=sigo(:,:,1);
polcoms.sig=sig;
polcoms.sigv=sigv;
polcoms.sigo=sigo;
polcoms.sigov=sigov;
polcoms.ds=ds;
polcoms.dsu=dsu;
polcoms.dsv=dsv;
polcoms.dsuv=dsuv;

return
