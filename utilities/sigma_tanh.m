function dist = sigma_tanh(nlev,dl,du)

  z=0.0;
  for k=1:nlev
    x1=dl+du;
    x1=x1*(nlev-1-k)/(nlev-1);
    x1=x1-dl;
    x1=tanh(x1);
    x2=tanh(dl);
    x3=x2+tanh(du);
    z(k+1)=(x1+x2)/x3-1.0;
  end;

