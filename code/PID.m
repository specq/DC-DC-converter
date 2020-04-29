clc;
clear;
close all;

% Electrical Parameters 
Rc = 0.33;
Rl = 1.6;
Ron = 0.05;
Ro = 100;
C = 56e-6;
L = 10e-3;
Vin = 15;
Vj = 0.1;

% Sampling period
Ts = 1e-3;

% PID parameters
Ki = 0.02;

% Computing the equilibrium point
x2_eq = 0;
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
C = [0 0.26];
y_eq = C*x_eq;

ss_c = ss(A,B,C,0);
ss_d = c2d(ss_c,Ts);

% Converting the system into transfer function
G_d = tf(ss_d);

% Definition of the discrete controller
z = tf('z',Ts);
K_d = Ki/(1-z^(-1));

% Creating the feedback loop
sys_fb = feedback(K_d*G_d,1);

% Poles and zeros
pzmap(sys_fb,'k');

% Step response
figure;
step(sys_fb);