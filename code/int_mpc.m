clc;
clear;
close all;

load ohm100_10kHz.mat;
G1 = ohm100_10kHz;
Ts = G1.Ts;

%A = (G1.A+G2.A)/2; B = (G1.B+G2.B)/2; C = [0 1];
A = G1.A;
A(1,2) = A(1,2)*100;
A(2,1) = A(2,1)/100;
B = G1.B; 
B(1) = B(1)*100;
C = [0 1];
ref = 5000;

Q = [1,0;0,10];
R = 1;

val_ss = [A-eye(2) B; C 0]\[0;0;1]*ref;
xs = val_ss(1:2);
us = round(val_ss(3));

% Horizon
N = 10;

% Defining the model
model = LTISystem('A',A,'B',B,'Ts',Ts);

%State constraints
model.x.min = [-xs(1); -xs(2)];
model.x.max = [20000-xs(1); 10000-xs(2)];

% Input constraints
model.u.min = -us;
model.u.max = 1000-us;

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

%expmpc = FittingController(expmpc);
expmpc.optimizer.trimFunction('primal', 1);
expmpc.optimizer.toMatlab('exp_sol_int.m', 'primal', 'first-region');

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

K = dlqr(A,B,[1,0;0,1],1);

x_hist = zeros(2,length(t));
u_hist = zeros(1,length(t)-1);
regions = zeros(1,length(t)-1);

for i=1:length(t)-1
    u = [];
    if x_hist(1,i+1)<21000 && x_hist(1,i+1)>20000
        [u,regions(i)] = exp_sol_int([20000;x_hist(2,i)]-xs);
    else
        [u,regions(i)] = exp_sol_int(x_hist(:,i)-xs);
    end
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
    if u_hist(i) > 1000
        u_hist(i) = 1000;
    end
    
    x_hist(:,i+1) = fix(A*x_hist(:,i)+B*u_hist(i));
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