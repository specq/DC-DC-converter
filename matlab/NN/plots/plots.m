clc;
clear;
close all;

load ohm100_10kHz.mat;
load param.mat
load controller.mat
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

t = 0:0.1:16;

K = dlqr(A,B,[1000000,0;0,1],1000);

x_sim = zeros(2,length(t));
x_sim_sc = zeros(2,length(t)-1);
u_sim = zeros(1,length(t)-1);
reg_sim = zeros(1,length(t)-1);

for i=1:length(t)-1
    dx = x_sim(:,i)-xs;
    x_sim_sc(:,i) = (dx - xmin) ./ xdelta;
    [u_temp,reg_sim(i)] = pwa(x_sim_sc(:,i));
    if isnan(u_temp)
        u_sim(i) = -K*dx;
        reg_sim(i)=0;
    else
        u = fit.O*u_temp + fit.o;
        u = (u*udelta) + umin;
        u_sim(i) = max(min(u,umax),umin);
    end
    u_sim(i) = u_sim(i) + us;
    
    x_sim(:,i+1) = A*x_sim(:,i)+B*u_sim(i);
end
x_sim(:,end)=[];
t(end) = [];

%% Experimental results
x1_exp = readmatrix('x1.csv');
x2_exp = readmatrix('x2.csv');
u_exp = readmatrix('u.csv');
reg_exp = readmatrix('regions.csv');

x1_exp(x1_exp(:,1)<0,:) = [];
x2_exp(x2_exp(:,1)<0,:) = [];
u_exp(u_exp(:,1)<0,:) = [];
reg_exp(reg_exp(:,1)<0,:) = [];

% Sampling of the scope data
x1_d = x1_exp(1,2);
x2_d = x2_exp(1,2);
u_d = [];
reg_d = [];

% Time of the scope data
time = 1000*u_exp(:,1);

% Sampling t for the state
T1 = 0.1;
% Sampling t for the duty cycle and regions
T2 = 0.05;
% Sampling
for i=1:length(time)
    if time(i)>T1
        x1_d = [x1_d;x1_exp(i,2)];
        x2_d = [x2_d;x2_exp(i,2)];
        T1 = T1 + 0.1;
    end
    if time(i)>T2
        u_d = [u_d;u_exp(i,2)];
        reg_d = [reg_d;reg_exp(i,2)];
        T2 = T2 + 0.1;
    end
end

% Shrinking and processing the sampled data
u_d=u_d(1:160,:);

% Computing the regions
reg_d=reg_d(1:160,:);
reg_exp = zeros(160,1);
reg_exp(reg_d<0.25) = 0; 
reg_exp(reg_d>0.25 & reg_d<0.75) = 1;
reg_exp(reg_d>0.75 & reg_d<1.25) = 2;
reg_exp(reg_d>1.25 & reg_d<1.75) = 3;
reg_exp(reg_d>1.75 & reg_d<2.25) = 4;
reg_exp(reg_d>2.25 & reg_d<2.75) = 5;

% Computing x1
x1_d = x1_d(1:160,:); 
x1_d = x1_d/7.5;

% Shrinking x2 table 
x2_d=x2_d(1:160,:);

% Scaling the state value 
x1_scaled = x1_d/xdelta(1);
x2_scaled = x2_d/xdelta(2);

% State plots
figure;
plot(mptSol.Set);
hold on;
p1=plot(x1_scaled, x2_scaled, 'k--o');
p2=plot(x_sim_sc(1,:),x_sim_sc(2,:),'b-o');
legend([p1;p2],{'Real system';'Simulation'});
ylabel('x_2 normalized');
xlabel('x_1 normalized');


% Time plots
figure;
subplot(4,1,1);
plot(t,1000*x1_d);
hold on;
plot(t,1000*x_sim(1,:));
grid on;
ylabel('x_1(mA)');

subplot(4,1,2);
plot(t,x2_d);
hold on;
plot(t,x_sim(2,:));
grid on;
ylabel('x_2(V)');

subplot(4,1,3);
plot(t,100*u_d);
hold on;
plot(t,100*u_sim);
grid on;
ylabel('Duty cycle(%)');
legend('Real system','Simulation');

subplot(4,1,4);
plot(t,reg_exp);
hold on;
plot(t,reg_sim);
grid on;
ylabel('Region ID');
xlabel('Time(ms)');