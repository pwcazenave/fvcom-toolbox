function dist = sigma_geo(nlev,p_sigma)

kb = nlev;
if(p_sigma ==1)
  for k=1:nlev
    dist(k) = -((k-1)/float(kb-1))^p_sigma;
  end;
else
  for k=1:(kb+1)/2
    dist(k) = -((k-1)/float((kb+1)/2-1))^p_sigma/2;
  end;
  do k=(kb+1)/2+1:kb   
    dist(k) = ((kb-k)/float((kb+1)/2-1))^p_sigma/2-1.0;
  end;
end;

