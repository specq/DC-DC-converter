clc;
clear;
close all;

% System dynamics
load ohm18
load ohm22
load ohm59
load ohm100
sys = ohm100;
A=sys.A;B=sys.B;C=sys.C;D=sys.D;
L = place(A',C',[0.1,0.11])'
Q = [100 0;0 1];
R = 1;
K = dlqr(A,B,Q,R);

ref = 5;

% Steady state parameters
N = [A-eye(2) B; C 0]\[0;0;1];
Nx = N([1,2]);
Nu = N(3);
offset = (Nu + K *Nx)*ref;
us = Nu*ref
xs = Nx*ref

%% Simulation
ref = 5;

t = 0:50;

x_hat = zeros(2,length(t));
x = zeros(2,length(t));
u = zeros(1,length(t)-1);
x_hat(:,1) = [1;2];
x(:,1) = [0;0];
integral = 0;

for i=1:length(t)-1
    integral = integral + ohm100.Ts*(ref-ohm100.C*x(:,i));
    u(i) = -K*x_hat(:,i) + offset + 1.5*integral;
    x_hat(:,i+1) = ohm100.A*x_hat(:,i)+ohm100.B*u(i)+L*ohm100.C*(x(:,i)-x_hat(:,i));
    x(:,i+1) = ohm100.A*x(:,i) + ohm100.B*u(i);
end

figure;
subplot(3,1,1);
grid on; hold on;
plot(t, x(1,:),'r');
plot(t, x_hat(1,:),'b');
grid on;

subplot(3,1,2);
grid on; hold on;
plot(t, x(2,:),'r');
plot(t, x_hat(2,:),'b');
legend('x','x_{hat}');
grid on;

subplot(3,1,3);
grid on; hold on;
plot(t(1:end-1), u);
grid on;

