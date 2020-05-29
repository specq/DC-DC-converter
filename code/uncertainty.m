clc;
clear;
%close all;

load ohm100_10kHz.mat;
load ohm50_10kHz.mat;
G1 = ohm100_10kHz;
G2 = ohm50_10kHz;
Ts = G1.Ts;

%A = (G1.A+G2.A)/2; B = (G1.B+G2.B)/2; C = [0 1];
A = G1.A; B = G1.B; C = [0 1];
ref = 5;

Q = [10000,0;0,1];
R = 1;

K = dlqr(A,B,Q,R);
val_ss = [A-eye(2) B; C 0]\[0;0;1]*ref;
xs = val_ss(1:2);
us = val_ss(3);
 
%% Simulation

t = 0:0.1:20;

x_hist = zeros(2,length(t));
u_hist = zeros(1,length(t)-1);
integral = 0;

for i=1:length(t)-1
    error = ref-C*x_hist(:,i);
    
    if error < -0.05 && integral > 0
        integral = 0;
    end
    
    if abs(error)<1
        integral = integral + error;
    else
        integral = 0;
    end
    
    u_hist(i) = -K*(x_hist(:,i)-xs)+us + 0.0*integral;
    
    if u_hist(i) < 0
        u_hist(i) = 0;
    end
    if u_hist(i) >1
        u_hist(i) = 1;
    end
    
    x_hist(:,i+1) = G1.A*x_hist(:,i)+G1.B*u_hist(i);
    %{
    if i <= 166
        x_hist(:,i+1) = G1.A*x_hist(:,i)+G1.B*u_hist(i);
    end
    if i > 166 && i <= 333
        x_hist(:,i+1) = G2.A*x_hist(:,i)+G2.B*u_hist(i);
    end
    if i > 333
        x_hist(:,i+1) = G1.A*x_hist(:,i)+G1.B*u_hist(i);
    end
    %}
    
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