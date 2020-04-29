clc;
clear;
close all;

load ohm100
nom = ohm100;
Ts = nom.Ts;
A = nom.A; B=nom.B; C = nom.C;
nx = size(B,1);
nu = size(B,2);
ny = size(C,1);

Q = [100 0; 0 1];
R = 1*eye(nu);

K = lqi(nom, [Q [0;0];0 0 1], R);
K0 = K([1,2]);
K1 = K(3);


N = [nom.A-eye(nx) nom.B; nom.C 0]\[0;0;1];
Nx = N([1,2]);
Nu = N(3);
N = Nu + K0 *Nx;
%% Simulation


ref = 5;
x0 = [0;0];
integral = 0;

t = 0:50;

x_hist = zeros(nx,length(t));
x_hist(:,1) = x0;

u_hist = zeros(1,length(t)-1);


for i=1:length(t)-1
    integral = integral + ohm100.Ts*(ref-C*x_hist(:,i));
    u_hist(:,i) = -K0*x_hist(:,i) + N*ref + K1 * integral;
    x_hist(:,i+1) = A * x_hist(:,i) + B * u_hist(:,i);
end


figure;
subplot(3,1,1);
grid on; hold on;
plot(t, x_hist(1,:));
grid on;

subplot(3,1,2);
grid on; hold on;
plot(t, x_hist(2,:));
grid on;

subplot(3,1,3);
grid on; hold on;
plot(t, [0 u_hist]);
grid on;




