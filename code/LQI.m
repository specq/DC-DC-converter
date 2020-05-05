clc;
clear;
close all;

load ohm100
load ohm72
load ohm36
load ohm18
G1 = ohm100;
G2 = ohm72;
G3 = ohm36;
G4 = ohm18;
Ts = 0.001;

G = stack(1,G1,G2,G3,G4);
C = [0,1];

ref = 5;

Q(:,:,1) = [2000 0; 0 1];
Q(:,:,2) = [2000 0; 0 1];
Q(:,:,3) = [10 0; 0 1];
Q(:,:,4) = [10 0; 0 1];
R = 1;

for i=1:4
    K(i,:) = dlqr(G(:,:,i,1).A,G(:,:,i,1).B,Q(:,:,i),R);
    val_ss(:,i) = [G(:,:,i,1).A-eye(2) G(:,:,i,1).B; C 0]\[0;0;1]*ref;
    xs(:,i) = val_ss(1:2,i);
    us(:,i) = val_ss(3,i);
end  
%% Simulation

Tf = 101;
t = 0:Tf-1;

x_hist = zeros(2,Tf);
x_hist(:,1) = [0;0];
u_hist = zeros(1,Tf-1);

% Actual model
A = ohm100.A; B = ohm100.B;

% First input
u_hist(1) = -K(1,:)*(x_hist(:,1)-xs(:,1))+us(1);

for i=2:Tf-1
    % Real system
    x_hist(:,i) = A*x_hist(:,i-1)+B*u_hist(i-1);
    y = C*x_hist(:,i);
    
    % Estimates
    error = [];
    for j=1:4
        error(j) = abs(y-C*(G(:,:,j,1).A*x_hist(:,i-1)+...
                           G(:,:,j,1).B*u_hist(i-1)));
    end
    % Chose the best model
    [~,sigma] = min(error);
    
    % Apply the input
    u_hist(i) = -K(sigma,:)*(x_hist(:,i)-xs(:,sigma)) + us(sigma);
    
    % Load switch
    if i == 33
        A = ohm18.A; B = ohm18.B;
    end
    if i == 66
        A = ohm100.A; B = ohm100.B;
    end
end
x_hist(:,Tf) = A*x_hist(:,Tf-1)+B*u_hist(Tf-1);


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