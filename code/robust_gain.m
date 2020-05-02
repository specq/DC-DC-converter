clc
clear
close all

load ohm100_c
load ohm22_c
A1 = ohm100_c.A; B1 = ohm100_c.B;
A2 = ohm22_c.A; B2 = ohm22_c.B;

nx = 2;
nu = 1;

Y = sdpvar(nu,nx);
L = sdpvar(nx,nx,'symmetric');

constraints = [L*A1'+A1*L-B1*Y-Y'*B1'+B1*B1' <= 0,...
               L*A2'+A2*L-B2*Y-Y'*B2'+B2*B2' <= 0,...
               L >= 0];
opt = sdpsettings('solver','mosek');
optimize(constraints,trace(L),opt);
L = value(L);
Y = value(Y);
K=Y/L;
sys1 = ss(A1-B1*K, B1, eye(2),[0;0]);
sys2 = ss(A2-B2*K, B2, eye(2),[0;0]);
G = stack(1,sys1,sys2);
step(G)
