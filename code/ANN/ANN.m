clc
close all
clear

load data_100ohms

X = data_100ohms(1:2,:);
Y = data_100ohms(3,:);

nz = 7;
nx = 2;
F = zeros(nz,nx);
f = zeros(nz,1);
G = zeros(1,nz);
g = 0;
L = zeros(nz,nz);
eps = 0;





