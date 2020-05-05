clear
clc
close all

load ohm100.mat
sys = ohm100;
A=sys.A; B=sys.B; C = sys.C; D=sys.D;
ref = 5;

% Steady state target
ss = [A-eye(2) B;C 0]\[0;0;1]*ref;
xss = ss(1:2);
uss = ss(3);

model = LTISystem('A',A,'B',B,'Ts',6.6667e-4);

%State constraints
model.x.min = [-xss(1); -xss(2)];
model.x.max = [0.5-xss(1); 10-xss(2)];

% Input constraints
model.u.min = -uss;
model.u.max = 1-uss;

% State penalty
Q = [100 0; 0 1];
model.x.penalty = QuadFunction(Q);

% Input penalty
R = 1;
model.u.penalty = QuadFunction(R);

% Terminal set and cost
Tset = model.LQRSet;
PN = model.LQRPenalty;
model.x.with('terminalSet');
model.x.terminalSet = Tset;
model.x.with('terminalPenalty');
model.x.terminalPenalty = PN;

% Horizon
N = 10;
mpc = MPCController(model, N);


%% Explicit MPC
close all
expmpc = mpc.toExplicit();
figure;
plot(expmpc.optimizer.Set);
X1=-xss(1):0.01:0.5-xss(1);
X2=-xss(2):0.1:10-xss(2);
Z= zeros(length(X1),length(X2));
for i=1:length(X1)
    for j=1:length(X2)
        Z(i,j) = exp_sol([X1(i),X2(i)]);
    end
end
[X1,X2]=meshgrid(X1,X2);
surf(X1',X2',Z);
%expmpc.exportToC('eMPC','directory');
expmpc.optimizer.trimFunction('primal', 1);
expmpc.optimizer.toMatlab('exp_sol.m', 'primal', 'obj')

%% Simulation
x0 = -xss;

t = 0:20;
x_hist = zeros(2,length(t)); x_hist(:,1)=x0;
u_hist = zeros(1,length(t)-1);

for i = 1:length(t)-1
    uopt = exp_sol(x_hist(:,i));
    u_hist(:,i) = uopt(1);
    x_hist(:,i+1) = A*x_hist(:,i) + B*u_hist(:,i);
end

figure;
subplot(3,1,1);
plot(t,xss(1)+x_hist(1,:));
grid on;
ylabel('x_1(mA)');

subplot(3,1,2)
plot(t,xss(2)+x_hist(2,:));
grid on;
ylabel('x_2(V)');

subplot(3,1,3)
plot(t(1:end-1),uss+u_hist);
grid on;
ylabel('Duty cycle(%)');
xlabel('Time(ms)');

