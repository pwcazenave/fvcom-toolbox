function [metvar,X_send,Y_send] = extract_mesoscale(fname,ndays)

%load in the mesoscale data from mesoscale operational POLCOMS met input data and convert to cs3 operational surge model met data grid
% Script from JMB

A= load (fullfile('/work/kthurs/Met_data_processing_scripts/MESOSCALE/OUTPUT/',fname));
tint=24/3;%Time interval 3hrs
days=ndays+2-1/tint; % number of days in month +2-1 to include extra day of output up to 21:00, the output before midnight of the next day
               % day 1 = 00:00 of the start day, day 1.5 = 12:00 of the start day

%POLCOMS meso met grid
dx=0.11; dy=0.11; 
x1=-13;
y1=48.39;
nx=218;
ny=136;
x2=(nx-1)*dx+x1;
y2=(ny-1)*dy+y1;
[X,Y]=meshgrid(x1:dx:x2,y1:dy:y2);%grid
X_send = x1:dx:x2;
Y_send = y1:dy:y2;

C=zeros(ny*days*tint,nx);%ygrid*number of days * outputs per day, xgrid
jj=1;
for ii=1:(ny*nx):length(A) %total length every 3 hours output  
   % B = reshape(A(ii,ii+29647),136,218);
    C(jj:jj+ny-1,:) = reshape(A(ii:ii+nx*ny-1),nx,ny)';
    jj=jj+ny;
end

% Reshape C to have form ny x nx x days*tint
metvar = zeros(ny,nx,days*tint);
jj=1;
for kk = 1:ny:length(C)
    metvar(:,:,jj) = C(kk:kk+ny-1,:);
    jj = jj+1;
end

%%
% % Load coastline
%  load '~/POLCOMS/POLCOMS_matlab/N_Atlantic_coast_m.dat';
%  x=N_Atlantic_coast_m(:,1);
%  y=N_Atlantic_coast_m(:,2);
% 
% % 
% for kk=1:ny:length(C)
%     figure(3)
%     pcolor(X,Y,C(kk:kk+ny-1,:))
%     shading interp
%     hold on
%     colorbar
% %     %Plot coastline
% %     plot(x,y,'k')
% %     drawnow
% end


