
% convert admesh edge matlab file from cartesian to spherical
base_dir = '~/Models/git/fvcom/locate-implementation/loc_v0/grids/lymebay/admesh/'
addpath '/users/modellers/rito/Code/fvcom-toolbox/utilities/'
admesh_file = 'Locate_500m'

Mobj=read_admesh_mesh('msh',fullfile(base_dir,[admesh_file,'.14']),'coordinate','spherical');
write_SMS_2dm(fullfile(base_dir,[admesh_file,'.2dm']),Mobj.tri,Mobj.lon,Mobj.lat,Mobj.h,[])
    z1 = utmzone(50,-4);
    [ellipsoid,estr] = utmgeoid(z1);
    utmstruct = defaultm('utm');
    utmstruct.zone = z1;
    utmstruct.geoid = ellipsoid;
    utmstruct = defaultm(utmstruct);
%coast = load([base_dir 'coast_rosa.mat']);
[Mobj.x,Mobj.y]=mfwdtran(utmstruct,Mobj.lat,Mobj.lon);
Mobj.have_xy = 1;
Mobj.have_lonlat = 0;

[Mobj] = write_admesh_mesh(Mobj,'output_directory',base_dir,'native_coord','cartesian');

write_SMS_2dm(fullfile(base_dir,[admesh_file,'admesh.2dm']),Mobj.tri,Mobj.x,Mobj.y,Mobj.h,[])


% read cartesian 2dm mesh from sms to add y-shift appearing after admesh
% processing
Mobj = read_sms_mesh(...
    '2dm', fullfile('~/Models/git/fvcom/locate-implementation/loc_v0/grids/', ['lyme_bay_v05_post_admesh', '.2dm']),...
    'coordinate', 'cartesian', 'addCoriolis', false);
%     'bath', fullfile(conf.base,conf.casename, [conf.casename, '.dat']),...
