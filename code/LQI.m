clc;
clear;
close all;

load ohm100_1kHz.mat;
load ohm100_10kHz.mat;
load ohm18_1kHz.mat;
load ohm18_10kHz.mat;
G1 = ohm100_10kHz;
G2 = ohm18_10kHz;
%G3 = ohm36;
%G4 = ohm18;
Ts = G1.Ts;

G = stack(1,G1,G2);
C = [0,1];

ref = 5;

Q(:,:,1) = [100 0; 0 1];
Q(:,:,2) = [100 0; 0 1];
%Q(:,:,3) = [10 0; 0 1];
%Q(:,:,4) = [10 0; 0 1];
R = 20;

for i=1:2
    K(i,:) = dlqr(G(:,:,i,1).A,G(:,:,i,1).B,Q(:,:,i),R);
    val_ss(:,i) = [G(:,:,i,1).A-eye(2) G(:,:,i,1).B; C 0]\[0;0;1]*ref;
    xs(:,i) = val_ss(1:2,i);
    us(:,i) = val_ss(3,i);
end  
%% Simulation

t = 0:0.1:20;

x_hist = zeros(2,length(t));
u_hist = zeros(1,length(t)-1);

% Actual model
A = G2.A; B = G2.B;

% First input
u_hist(1) = -K(1,:)*(x_hist(:,1)-xs(:,1))+us(1);

for i=2:length(t)-1
    % Real system
    x_hist(:,i) = A*x_hist(:,i-1)+B*u_hist(i-1);
    y = C*x_hist(:,i);
    
    % Estimates
    error = [];
    for j=1:2
        error(j) = abs(y-C*(G(:,:,j,1).A*x_hist(:,i-1)+...
                            G(:,:,j,1).B*u_hist(i-1)));
    end
    % Chose the best model
    [~,sigma] = min(error);
    
    % Apply the input
    u_hist(i) = -K(sigma,:)*(x_hist(:,i)-xs(:,sigma)) + us(sigma);
    %{
    % Load switch
    if i == 33
        A = G2.A; B = G2.B;
    end
    if i == 66
        A = G1.A; B = G1.B;
    end
    %}
end
x_hist(:,length(t)) = A*x_hist(:,length(t)-1)+B*u_hist(length(t)-1);


figure;
subplot(3,1,1);
grid on; hold on;
plot(t, x_hist(1,:));
grid on;
ylabel('x_1(mA)');

subplot(3,1,2);
grid on; hold on;
plot(t, x_hist(2,:));
grid on;
ylabel('x_2(V)');

subplot(3,1,3);
grid on; hold on;
plot(t, [0 u_hist]);
ylabel('Duty cycle(%)');
xlabel('Time(ms)');
grid on;