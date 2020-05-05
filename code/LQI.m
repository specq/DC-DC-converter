clc;
clear;
close all;

load ohm100
load ohm18
sys = ohm100;
Ts = sys.Ts;
A = (ohm100.A+ohm18.A)./2; B=(ohm100.B+ohm18.B)./2; C = sys.C;
ref = 5;

Q = [100 0; 0 1];
R = 1;

K = dlqr(A,B,Q,R)

val_ss = [A-eye(2) B; C 0]\[0;0;1]*ref;
xs = val_ss(1:2)
us = val_ss(3)
%% Simulation

integral = 0;
Ki = 10;
Kp = 0.0;
t = 0:100;

x_hist = zeros(2,length(t));
x_hist(:,1) = [0;0];

u_hist = zeros(1,length(t)-1);
%A = ohm100.A;B=ohm100.B;

for i=1:length(t)-1
    error = ref-C*x_hist(:,i);
    integral = integral + Ts*error;
    u_hist(:,i) = -K*(x_hist(:,i)-xs) + us + Ki*integral + Kp*error;
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




