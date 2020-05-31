clc;
clear;
close all;

load ohm50_10kHz.mat;
load ohm100_10kHz.mat;
G1 = ohm100_10kHz;
G2 = ohm50_10kHz;
Ts = G1.Ts;

%A = (G1.A+G2.A)/2; B = (G1.B+G2.B)/2; C = [0 1];
A = G1.A; B = G1.B; C = [0 1];
ref = 5;

Q = [90,0;0,1];
R = 1;

val_ss = [A-eye(2) B; C 0]\[0;0;1]*ref;
xs = val_ss(1:2);
us = val_ss(3);

% Horizon
N = 10;

% Defining the model
model = LTISystem('A',A,'B',B,'Ts',Ts);

%State constraints
model.x.min = [-xs(1); -xs(2)];
model.x.max = [0.2-xs(1); 11.7-xs(2)];

% Input constraints
model.u.min = -us;
model.u.max = 1-us;

% State penalty
model.x.penalty = QuadFunction(Q);

% Input penalty
model.u.penalty = QuadFunction(R);


% Terminal set and cost
Tset = model.LQRSet;
PN = model.LQRPenalty;
model.x.with('terminalSet');
model.x.terminalSet = Tset;
model.x.with('terminalPenalty');
model.x.terminalPenalty = PN;

mpc = MPCController(model, N);

% Explicit MPC
expmpc = mpc.toExplicit();

expmpc = FittingController(expmpc);
expmpc.optimizer.trimFunction('primal', 1);
expmpc.optimizer.toMatlab('exp_sol.m', 'primal', 'first-region');

%% Plots 
figure;
expmpc.feedback.fplot()
H = {}; h = {};
F = {}; g = {};
for i = 1:expmpc.optimizer.Num
        F{i} = expmpc.optimizer.Set(i).Functions('primal').F;
        g{i} = expmpc.optimizer.Set(i).Functions('primal').g;
        H{i} = expmpc.optimizer.Set(i).A;
        h{i} = expmpc.optimizer.Set(i).b;
end
%% Simulation

t = 0:0.1:20;

K = dlqr(G1.A,G1.B,[1000,0;0,1],1000);

x_hist = zeros(2,length(t));
u_hist = zeros(1,length(t)-1);
regions = zeros(1,length(t)-1);
integral = 0;

for i=1:length(t)-1
    error = ref-C*x_hist(:,i);
   
    if error < -0.05 && integral > 0
        integral = 0;
    end
    
    integral = integral + error;
    
    [u,regions(i)] = exp_sol(x_hist(:,i)-xs);
    if isnan(u)
        u_hist(i) = -K*(x_hist(:,i)-xs);
        regions(i)=0;
    else
        u_hist(i) = u;
    end
    u_hist(i) = u_hist(i) + us;
    
    if u_hist(i) < 0
        u_hist(i) = 0;
    end
    if u_hist(i) >1
        u_hist(i) = 1;
    end
    
    x_hist(:,i+1) = G1.A*x_hist(:,i)+G1.B*u_hist(i);
end

figure;
plot(expmpc.optimizer.Set);
hold on;
plot(x_hist(1,:)-xs(1)*ones(1,length(t)), x_hist(2,:)-xs(2)*ones(1,length(t)), 'k--o');

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

%% 
A2 = round(1000*A);
B2 = round(1000*B);
val_ss = [A2-1000*eye(2) B2; C 0]\[0;0;1]*5000;
xs2 = [50;5000];
us2 = 351;
fix((A2*xs2+B2*us2)/1000)
