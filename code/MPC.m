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
x_eq = [x1_eq;5];
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
H = [1 0;-1 0;0 1];
h = [1; 1; 10];
G = [1;-1];
g = [1;0];

ss_c = ss(A,B,C,0);
ss_d = c2d(ss_c,Ts);

nx = size(B,1);
nu = size(B,2);
Q = eye(nx);
R = eye(nu);

opt = stepDataOptions('StepAmplitude',u_eq,'InputOffset',-u_eq);
[y,t]=step(ss_c,0.025,opt);
figure;
plot(1000*t,y+x2_eq);
xlabel('Time (ms)');
ylabel('x2 (V)');
grid on;
%% Terminal set
K = dlqr(ss_d.A,ss_d.B,Q,R);
K = -K;
Ak = ss_d.A + ss_d.B*K;

O = polytope([H;G*K],[h;g]);
figure;
p1 = plot(O);


while 1
    Oprev = O;
    [F,f]=double(O);
    O = polytope([F;F*Ak],[f;f]);
    if O == Oprev
        break;
    end
end
