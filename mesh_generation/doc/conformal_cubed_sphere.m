clc
clear

iPatch = 2;

ccs_mesh = '..\run\ccs_output.nc';
R2D = 180/pi;
D2R = pi/180;
x   = ncread(ccs_mesh,'x');
y   = ncread(ccs_mesh,'y');
lon = ncread(ccs_mesh,'lon');
lat = ncread(ccs_mesh,'lat');

lon1d = reshape(lon,[],1);
lat1d = reshape(lat,[],1);

radius = 1;
[xp,yp,zp] = sph2cart(lon1d*D2R,lat1d*D2R,radius);

% Plot white base
[xs,ys,zs] = sphere(800);
s = surf(xs*radius,ys*radius,zs*radius,'FaceColor','w');
set(s,'EdgeColor','None')
hold on

scatter3(xp,yp,zp,'.')

% figure
% plt = pcolor(squeeze(lat(:,:,iPatch)));
% set(plt,'EdgeColor','none')
% colormap(jet)

% patch(pointCoord(1,pointId_patch),pointCoord(2,pointId_patch),pointCoord(3,pointId_patch),areaCell(iCell));

% figure
% lon1d = reshape(lon(:,:,iPatch),[],1);
% lat1d = reshape(lat(:,:,iPatch),[],1);
% 
% [xp,yp,zp] = sph2cart(lon1d*D2R,lat1d*D2R,radius);
% scatter3(xp,yp,zp,'.')