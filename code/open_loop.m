clc;
clear;
close all;

% Electrical Parameters 
Rc = 0.33;
Rl = 2;
Ron = 0.05;
Ro = 100;
C = 56e-6;
L = 10e-3;
Vin = 15;
Vj = 0.1;

% Sampling period
Ts = 2e-4;

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
C = [1 0];

if det([C;C*A]) ~= 0 && det([B A*B]) ~= 0
    disp("Observable and controllable");
end

nx = size(B,1);
nu = size(B,2);
ny = size(C,1);

s=tf('s');
omega = 2*pi*100;
W = omega/(omega+s);

G = ss(A,B,C,0);
ohm100_5kHz = c2d(G,Ts);
