%! this program adds sigma levels/layers and shifted layers to a netcdf file 
%! so that Datatank can visualize the results.
%! 

function add_siglay_forDT(fname)

if(~exist(fname))
  error(['file: ' fname ' does not exist'])
end;
fprintf('adding siglay DT stuff to %s\n',fname);
nc = netcdf(fname,'w');

siglay = nc{'siglay'}(:,:);
siglev = nc{'siglev'}(:,:);
[nlay,node] = size(siglay);
[nlev,node] = size(siglev);
fprintf('siglay dimension is %d\n',nlay);
fprintf('siglev dimension is %d\n',nlev);

siglayDT = zeros(nlay,1);
siglay2DT = zeros(nlay,1);
siglevDT = zeros(nlev,1);

for i=1:nlev
  siglevDT(i) = double(i-1)/double(nlev-1);
end;
for i=1:nlay
  siglayDT(i) = 0.5*(siglevDT(i)+siglevDT(i+1));
  siglay2DT(i) = siglevDT(i); 
end;

nc{'siglay_DT'} = ncfloat('siglay') ;
nc{'siglay_DT'}.long_name = 'siglay_DT'; 

nc{'siglay_shift_DT'} = ncfloat('siglay') ;
nc{'siglay_shift_DT'}.long_name = 'siglay_shift_DT'; 

nc{'siglev_DT'} = ncfloat('siglev') ;
nc{'siglev_DT'}.long_name = 'siglev_DT'; 

nc{'siglay_DT'}(:) = siglayDT;
nc{'siglay_shift_DT'}(:) = siglay2DT;
nc{'siglev_DT'}(:) = siglevDT;

nc = close(nc);


