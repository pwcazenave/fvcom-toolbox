function dist = sigma_gen(nlev,dl,du,kl,ku,zkl,zku,h,hmin)


 
  if(h < hmin)
    z(1) = 0.0;
    dl2=0.001;
    du2=0.001;
    for k=1:nlev-1
      x1 = dl2+du2;
      x1 = x1*double(nlev-1-k)/double(nlev-1);
      x1 = x1 - dl2;
      x1 = tanh(x1);
      x2 = tanh(dl2);
      x3 = x2+tanh(du2);
      z(k+1) = (x1+x2)/x3-1.0;  
    end;
  else;
    %dr=(h-sum(zku)-sum(zkl))/h/double(nlev-ku-kl-1);
    dr=(h-du-dl)/h/double(nlev-ku-kl-1);
    z(1) = 0.0;
 
    for k=2:ku+1 
       z(k) = z(k-1)-zku(k-1)/h;
  %     fprintf('building z %f %f %f %f \n',z(k),zku(k-1),h,zku(k-1)/h)
    end;

    for k=ku+2:nlev-kl
      z(k)=z(k-1)-dr;
  %     fprintf('building z %f %f \n',z(k),dr)
    end;

    kk = 0;
    for k=nlev-kl+1:nlev
      kk=kk+1;
      z(k)=z(k-1)-zkl(kk)/h;
  %     fprintf('building z %f %f \n',z(k),zkl(kk))
    end;   
  end;
  dist = z;
