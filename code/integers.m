clc;
clear;
close all;

load ohm100_10kHz.mat;
G1 = ohm100_10kHz;
Ts = G1.Ts;

%A = (G1.A+G2.A)/2; B = (G1.B+G2.B)/2; C = [0 1];
A = G1.A; B = G1.B; C = [0 1];
A(1,2) = A(1,2)*100;
A(2,1) = A(2,1)/100;
B(1) = B(1)*100;
ref = 5000;

Q = [10000,0;0,1];
R = 1;

K = dlqr(A,B,Q,R);
val_ss = [A-eye(2) B; C 0]\[0;0;1]*ref;
xs = val_ss(1:2);
us = round(val_ss(3));
 
%% Simulation

t = 0:0.1:50;

x_hist = zeros(2,length(t));
u_hist = zeros(1,length(t)-1);

for i=1:length(t)-1
    u_hist(i) = us;
    x_hist(:,i+1) = fix(A*x_hist(:,i)+B*u_hist(i)); 
    if x_hist(1,i+1) < 0 
        x_hist(1,i+1) = 0;
    end
end


figure;
subplot(3,1,1);
grid on; hold on;
plot(t, x_hist(1,:)/100);
ylabel('x_1(mA)');

subplot(3,1,2);
grid on; hold on;
plot(t, x_hist(2,:)/1000);
ylabel('x_2(V)');

subplot(3,1,3);
grid on; hold on;
plot(t, [0 u_hist]/10);
ylabel('Duty cycle(%)');
xlabel('Time(ms)');