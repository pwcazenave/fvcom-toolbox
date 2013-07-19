function [pts] = inkmlbb(fname,x,y)
%fname = 'Wet_Skagit_BB.kml';
%x = geog(:,1);
%y = geog(:,2);


triangle = '/usr/local/bin/triangle';

[latb,lonb,dumz] = read_kml(fname);
[ktri,kx,ky] = tri_from_bndry(lonb,latb,triangle);

box = [min(kx),max(kx),min(ky),max(ky)];
pts1 = inbox(box,x,y);

mark = zeros(length(x),1);

dims = size(ktri);

for ii=1:length(pts1);
i = pts1(ii);
for j=1:dims(1);
  xtri = kx(ktri(j,1:3));
  ytri = ky(ktri(j,1:3));
  if(isintriangle(xtri,ytri,x(i),y(i)) ) ; mark(i) = 1; end;
end;
end;

pts = find(mark==1);

%figure
%patch('Vertices',[kx(:),ky(:)],'Faces',ktri,...
%       'Cdata',kx,'edgecolor','k','facecolor','k');
%hold on;
%plot(x(pts1),y(pts1),'r+');
%plot(x(pts),y(pts),'g+');
%axis equal
%colorbar

