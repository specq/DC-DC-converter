clc;
clear;
close all;

% Electrical Parameters 
Rc = 0.33;
Rl = 1.6;
Ron = 0.05;
Ro = 22;
C = 56e-6;
L = 10e-3;
Vin = 15;
Vj = 0.1;

% Sampling period
Ts = 1e-3;

% Computing the equilibrium point
x2_eq = 5;
x1_eq = x2_eq/Ro;
x_eq = [x1_eq;x2_eq];
u_eq = (Ro*Vj + Rl*x2_eq + Ro*x2_eq)/(Ro*Vin + Ro*Vj - Ron*x2_eq);

% Creating the state space model
a11 = -(Rl+Ron*u_eq)/L;
a12 = -1/L;
a21 = (-Rc*Ro*Rl*C+Ro*L-Rc*Ro*Ron*C*u_eq)/((Rc+Ro)*L*C);
a22 = -(Rc*Ro*C+L)/((Rc+Ro)*L*C);
b1 = (Vin+Vj-Ron*x1_eq)/L;
b2 = (Rc*Ro*(Vin+Vj-Ron*x1_eq))/((Rc+Ro)*L);

A = [a11 a12; a21 a22];
B = [b1;b2];
C = [0 1];

if det([C;C*A]) ~= 0 && det([B A*B]) ~= 0
    disp("Observable and controllable");
end

nx = size(B,1);
nu = size(B,2);
ny = size(C,1);

ss_c = ss(A,B,C,0);
ss_d = c2d(ss_c,Ts);

Q = [100 0; 0 1];
R = 1*eye(nu);
K_continuous = lqr(ss_c.A,ss_c.B,Q,R);
K_continuous = -K_continuous
K_discrete = dlqr(ss_d.A,ss_d.B,Q,R);
K_discrete = -K_discrete;

%% Simulation

ref = 5;
x0 = [0;0];

t = 0:50;

x_hist = zeros(nx,length(t));
x_hist(:,1) = x0;

u_hist = zeros(1,length(t)-1);

for i=1:length(t)-1
    u_hist(:,i) = K_discrete*(x_hist(:,i)-x_eq) + u_eq;
    x_hist(:,i+1) = x_eq + ss_d.A * (x_hist(:,i)-x_eq) + ss_d.B * (u_hist(:,i)-u_eq);
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



