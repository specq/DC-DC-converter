clc;
clear;
close all;

load ohm100_10kHz.mat;
load param.mat
load simplified_mpc_nn.mat
G = ohm100_10kHz;
Ts = G.Ts;
ref = 5;

A = G.A; B = G.B; C = [0 1];

val_ss = [A-eye(2) B; C 0]\[0;0;1]*ref;
xs = val_ss(1:2);
us = val_ss(3);

xmax = [0.15; 2]; 
xmin = [-0.05; -5];

umax = 0.6621;
umin = -0.3379;

xdelta = xmax-xmin;
udelta = umax-umin;

%% Simulation

t = 0:0.1:20;

K = dlqr(A,B,[1000000,0;0,1],1000);

x_hist = zeros(2,length(t));
u_hist = zeros(1,length(t)-1);
regions = zeros(1,length(t)-1);

for i=1:length(t)-1
    dx = x_hist(:,i)-xs;
    temp = (dx - xmin) ./ xdelta;
    [u_temp,regions(i)] = fitMpc(temp);
    if isnan(u_temp)
        u_hist(i) = -K*dx;
        regions(i)=0;
    else
        u = fit.O*u_temp + fit.o;
        u = (u*udelta) + umin;
        u_hist(i) = max(min(u,umax),umin);
    end
    u_hist(i) = u_hist(i) + us;
    
    x_hist(:,i+1) = A*x_hist(:,i)+B*u_hist(i);
end

figure;
subplot(3,1,1);
grid on; hold on;
plot(t, x_hist(1,:));
ylabel('x_1(mA)');

subplot(3,1,2);
grid on; hold on;
plot(t, x_hist(2,:));
ylabel('x_2(V)');

subplot(3,1,3);
grid on; hold on;
plot(t, [0 u_hist]);
ylabel('Duty cycle(%)');
xlabel('Time(ms)');

%% Get the parameters

H = {}; h = {};
F = {}; g = {};
for i = 1:mptSol.Num
        F{i} = mptSol.Set(i).Functions('primal').F;
        g{i} = mptSol.Set(i).Functions('primal').g;
        H{i} = mptSol.Set(i).A;
        h{i} = mptSol.Set(i).b;
end