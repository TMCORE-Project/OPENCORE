clc
clear

iPatch = 2;

ccs_mesh = '..\run\ccs_output.nc';
R2D = 180/pi;
D2R = pi/180;
x   = ncread(ccs_mesh,'x');
y   = ncread(ccs_mesh,'y');
z   = ncread(ccs_mesh,'z');
lon = ncread(ccs_mesh,'lon');
lat = ncread(ccs_mesh,'lat');

var = lon;

x1 = squeeze(var(:,:,1));
x2 = squeeze(var(:,:,2));
x3 = squeeze(var(:,:,3));
x4 = squeeze(var(:,:,4));
x5 = squeeze(var(:,:,5));
x6 = squeeze(var(:,:,6));