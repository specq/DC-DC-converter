clc;
clear;
close all;

% System dynamics
load ohm18
load ohm100

sys = ohm100;
A=sys.A;B=sys.B;C=eye(2);D=[0;0];Ts=sys.Ts;

L = place(A',C',[0.2,0.21,0.22])';

% Steady state target
ref = 5;
ss = [A-eye(2) B;0 1 0]\[0;0;1]*ref;
xss = ss(1:2);
uss = ss(3);

model = LTISystem('A',A,'B',B,'Ts',Ts);

%State constraints
model.x.min = [-xss(1); -xss(2)];
model.x.max = [0.5-xss(1); 10-xss(2)];

% Input constraints
model.u.min = -uss;
model.u.max = 1-uss;

% State penalty
Q = [100,0;0,1];
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

%% Simulation

t = 0:50;
x_hist = zeros(2,length(t)); x_hist(:,1)=-xss;
d_hist = zeros(2,length(t)); d_hist(:,1)=[0;0];
u_hist = zeros(1,length(t)-1);

x = -xss;
for i = 1:length(t)-1
    uopt = mpc.evaluate(x_hist(:,i));
    u_hist(:,i) = uopt(1);
    x_aug = A_aug*[x_hist(:,i);d_hist(:,i)] + B_aug*u_hist(:,i)...
            + L*(x-x_hist(1:2,i));
    x_hist(:,i+1) = x_aug(1:2);
    d_hist(:,i+1) = x_aug(3:4);
    x = A*x+B*u_hist(:,i);
end

figure;
subplot(5,1,1);
plot(t,x_hist(1,:));
grid on;
ylabel('x_1(mA)');

subplot(5,1,2)
plot(t,x_hist(2,:));
grid on;
ylabel('x_2(V)');

subplot(5,1,3)
plot(t,d_hist(1,:));
grid on;
ylabel('d1');

subplot(5,1,4)
plot(t,d_hist(2,:));
grid on;
ylabel('d2');

subplot(5,1,5)
plot(t(1:end-1),u_hist);
grid on;
ylabel('Duty cycle(%)');
xlabel('Time(ms)');
