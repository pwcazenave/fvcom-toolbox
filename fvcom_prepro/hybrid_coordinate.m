% calculate and write hybrid coordinate layer
function Mobj=hybrid_coordinate(conf,Mobj)
% defaults
% conf.sigma_file='coord_hybrid.sig'
% conf.nsigma=41;DU=25;DL=25;Hmin=200;KU=5;KL=5;ZKU=[.5 .5 .5 .5 .5];ZKL=[.5 .5 .5 .5 .5 ];
% conf.H0=100;conf.nlev=20;conf.DU=25;conf.DL=25;conf.KU=5;conf.KL=5;
myset=optimset('MaxFunEvals',5000,'MaxIter',5000,'TolFun',10e-5,'TolX',1e-7);
% solve for z1-z2 to find Hmin parameter
H0=conf.H0;
nlev=conf.nlev;DU=conf.DU;DL=conf.DL;
KU=conf.KU;KL=conf.KL;
ZKU=repmat(DU./KU,1,KU);ZKL=repmat(DL./KL,1,KL);
myfunparams=@(H)hybrid_coordinate_hmin(H,nlev,DU,DL,KU,KL,ZKU,ZKL);
% myfunparams(200)
[Hmin,fval] = fminsearch(myfunparams,H0,myset);
% print hybrid_coord.sigma file
fout=fopen(conf.sigma_file,'wt');
fprintf(fout,'%s\n',['NUMBER OF SIGMA LEVELS = ',num2str(nlev)]);
fprintf(fout,'%s\n',['SIGMA COORDINATE TYPE = GENERALIZED']);
fprintf(fout,'DU = %4.1f\n',(DU));
fprintf(fout,'DL = %4.1f\n',(DL));
fprintf(fout,'MIN CONSTANT DEPTH = %10.1f\n',round(Hmin));
fprintf(fout,'%s\n',['KU = ',num2str(KU)]);
fprintf(fout,'%s\n',['KL = ',num2str(KL)]);
fprintf(fout,'ZKU = ');
for ii=1:length(ZKU)
    fprintf(fout,'%4.1f',ZKU(ii));
end
fprintf(fout,'\n');
fprintf(fout,'ZKL = ');
for ii=1:length(ZKU)
    fprintf(fout,'%4.1f',ZKL(ii));
end
fprintf(fout,'\n');
fclose(fout);

for xx=1:length(Mobj.h)
    % loop through all nodes to create sig_coor
 Mobj.siglev(xx,:) = sigma_gen(nlev,DL,DU,KL,KU,ZKL,ZKU,Mobj.h(xx),Hmin);
end
Mobj.siglay = Mobj.siglev(:,1:end-1) + (diff(Mobj.siglev,1,2)/2);
for zz=1:nlev-1
Mobj.siglevc(:,zz) = nodes2elems(Mobj.siglev(:,zz), Mobj);
% Create a siglay variable (i.e. midpoint in the sigma levels).
Mobj.siglayc(:,zz) = nodes2elems(Mobj.siglay(:,zz), Mobj);
end
Mobj.siglevc(:,zz+1) = nodes2elems(Mobj.siglev(:,zz+1), Mobj);

% Create a depth array for the element centres.
if ~isfield(Mobj,'hc')
Mobj.hc = nodes2elems(Mobj.h, Mobj);
end

Mobj.siglevz = repmat(Mobj.h, 1, nlev) .* Mobj.siglev;
Mobj.siglayz = repmat(Mobj.h, 1, nlev-1) .* Mobj.siglay;
Mobj.siglevzc = repmat(Mobj.hc, 1, nlev) .* Mobj.siglevc;
Mobj.siglayzc = repmat(Mobj.hc, 1, nlev-1) .* Mobj.siglayc;

% % Create a siglay variable (i.e. midpoint in the sigma levels).
% zlay = z(1:end-1) + (diff(z)/2);
% 
% % Create a depth array for the element centres.
% hc = nodes2elems(Mobj.h, Mobj);
% 
% Mobj.siglevz = repmat(Mobj.h, 1, nlev) .* repmat(z, Mobj.nVerts, 1);
% Mobj.siglayz = repmat(Mobj.h, 1, nlev-1) .* repmat(zlay, Mobj.nVerts, 1);
% Mobj.siglevzc = repmat(hc, 1, nlev) .* repmat(z, Mobj.nElems, 1);
% Mobj.siglayzc = repmat(hc, 1, nlev-1) .* repmat(zlay, Mobj.nElems, 1);
% 
% % Add the sigma levels and layers to the Mobj.
% Mobj.siglev = z;
% Mobj.siglay = zlay;



return

% 
% % test with made up bathymetr
% Hmax=200;
% y=0:100;B=100;
% H = Hmax .*exp( -((y./B-0.15).^2 ./(0.5).^2) );

% H=120
% Z0=zeros(length(H),nsigma);Z2=zeros(length(H),nsigma);
% for hh=1:length(H)
%     Z0(hh,1)=0;DL2=0.001;DU2=0.001;KBM1=nsigma-1;
%     for nn=1:nsigma-1
%         X1=DL2+DU2;
%         X1=X1*(KBM1-nn)/KBM1;
%         X1=X1-DL2;
%         X1=tanh(X1);
%         X2=tanh(DL2);
%         X3=X2+tanh(DU2);
%         
%         Z0(hh,nn+1)=((X1+X2)/X3)-1.0;
%     end
%     
%     % s-coordinates
%     X1=(H(hh)-DU-DL);
%     X2=X1./H(hh);
%     DR=X2./(nsigma-KU-KL-1);
%     
%     Z2(hh,1)=0.0;
%     
%     for K=2:KU+1
%         Z2(hh,K)=Z2(hh,K-1)-(ZKU(K-1)./H(hh));
%     end
%     
%     for K=KU+2:nsigma-KL
%         Z2(hh,K)=Z2(hh,K-1)-DR;
%     end
%     
%     KK=0;
%     for K=nsigma-KL+1:nsigma
%         KK=KK+1;
%         Z2(hh,K)=Z2(hh,K-1)-(ZKL(KK)./H(hh));
%     end
% end
