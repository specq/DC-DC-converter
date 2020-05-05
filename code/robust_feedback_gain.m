clc;
close all;
clear;

load ohm100_c.mat
load ohm22_c.mat
G1 = ohm100_c;
G2 = ohm22_c;
G1_d = c2d(G1,0.001);
G2_d = c2d(G2,0.001);

A1 = G1.A; B1 = G1.B; 
A2 = G2.A; B2 = G2.B;

G = stack(1,G1,G2);

% Decision variables
L = sdpvar(2,2,'symmetric');
Y = sdpvar(1,2);
gamma1 = sdpvar(1,1);
gamma2 = sdpvar(1,1);


% Objective function
obj = gamma1+gamma2;

% Constraints
lmi = A1*L+L*A1'-B1*Y-Y'*B1'+B1*B1' <= 0;
lmi = [lmi, A2*L+L*A2'-B2*Y-Y'*B2'+B2*B2' <= 0];
lmi = [lmi,gamma1*eye(2)-L >= 0];
lmi = [lmi,[gamma2,Y;Y',L] >=0];
lmi = [lmi,L>=0];

% Solve the problem
options = sdpsettings('solver','mosek');
optimize(lmi,obj,options);
L = value(L);
Y = value(Y);

% Feedback gain matrix
K = Y / L;
sys_c = ss(A1-B1*K,B1,eye(2),[0;0]);
sys_d = c2d(sys_c,0.001);
K_d = -pinv(G1_d.B)*(sys_d.A-G1_d.A)';
sys_d.A-(G1_d.A-G1_d.B*K_d)

%% Simulation
x0 = [0;0];
ref = 5;
vs = [A1-eye(2) B1;0 1 0]\[0;0;1]*ref;
xs = vs(1:2);
us = vs(3);

t = 0:50;
x_hist = zeros(2,length(t)); x_hist(:,1)=x0;
u_hist = zeros(1,length(t)-1);

for i = 1:length(t)-1
    x_hist(:,i+1) = sys_d.A+x_hist(:,i);
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


