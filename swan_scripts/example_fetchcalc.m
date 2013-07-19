clear all; 
close all;
% determine fetch from Raubenheimer & Elgar Instrument 41 
re41fetch = fetch_calc('skg4.3_grd.dat','skg4.3_dep.dat',-122.4721646,48.3372476,-3,3,30,16,0.1,125);

% plot fetch 
show_fetch(re41fetch);

% calculate fetch for a west wind with 0.5m of water at the site
uwind = 10;
vwind = 0;
fprintf('with 0.5 meters of water\n')
[fetch] = get_fetch(re41fetch,uwind,vwind,0.5);
fprintf('fetch from west wind %f\n',fetch);
[fetch] = get_fetch(re41fetch,0,10,0.5);
fprintf('fetch from south wind %f\n',fetch);
[fetch] = get_fetch(re41fetch,-1,0,0.5);
fprintf('fetch from east wind %f\n',fetch);
[fetch] = get_fetch(re41fetch,0,-1,0.5);
fprintf('fetch from north wind %f\n',fetch);
fprintf('with 2.0 meters of water\n')
[fetch] = get_fetch(re41fetch,uwind,vwind,2.0);
fprintf('fetch from west wind %f\n',fetch);
[fetch] = get_fetch(re41fetch,0,10,2.0);
fprintf('fetch from south wind %f\n',fetch);
[fetch] = get_fetch(re41fetch,-1,0,2.0);
fprintf('fetch from east wind %f\n',fetch);
[fetch] = get_fetch(re41fetch,0,-1,2.0);
fprintf('fetch from north wind %f\n',fetch);
