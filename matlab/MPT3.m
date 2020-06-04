clear
clc
close all

load ohm100_1kHz.mat;
load ohm100_10kHz.mat;
load ohm18_1kHz.mat;
load ohm18_10kHz.mat;
load ohm100_5kHz.mat;
load ohm18_5kHz.mat;

G1 = ohm100_1kHz;
G2 = ohm18_1kHz;
Ts = G1.Ts;

G = stack(1,G1,G2);
C = [0,1];

ref = 5;

Q(:,:,1) = [100 0; 0 1];
Q(:,:,2) = [100 0; 0 1];
R = 1;

% Horizon
N = 10;

for i=1:2
    val_ss(:,i) = [G(:,:,i,1).A-eye(2) G(:,:,i,1).B; C 0]\[0;0;1]*ref;
    xs(:,i) = val_ss(1:2,i);
    us(i) = val_ss(3,i);
    
    % Defining the model
    model(i) = LTISystem('A',G(:,:,i,1).A,'B',G(:,:,i,1).B,'Ts',Ts);
    
    %State constraints
    model(i).x.min = [-xs(1,i); -xs(2,i)];
    model(i).x.max = [0.4-xs(1,i); 10-xs(2,i)];
    
    % Input constraints
    model(i).u.min = -us(i);
    model(i).u.max = 1-us(i);
    
    % State penalty
    model(i).x.penalty = QuadFunction(Q(:,:,i));
    
    % Input penalty
    model(i).u.penalty = QuadFunction(R);
    
    % Terminal set and cost
    Tset = model(i).LQRSet;
    PN = model(i).LQRPenalty;
    model(i).x.with('terminalSet');
    model(i).x.terminalSet = Tset;
    model(i).x.with('terminalPenalty');
    model(i).x.terminalPenalty = PN;

    mpc(i) = MPCController(model(i), N);
    
    % Explicit MPC
    expmpc(i) = mpc(i).toExplicit();
    %figure;
    %plot(expmpc(i).optimizer.Set);
    expmpc(i).optimizer.trimFunction('primal', 1);
    expmpc(i).optimizer.toMatlab(strcat('exp_sol_',int2str(i),'_1k.m'), 'primal', 'first-region');

end  
%% Simulation
x0 = [0;0];

t = 0:1:100;
x_hist = zeros(2,length(t)); 
u_hist = zeros(1,length(t)-1);

% Actual model
A = G1.A; B = G1.B;

% First input
u_hist(1) = exp_sol_1_1k(x_hist(:,1)-xs(:,1))+us(1);

for i = 2:length(t)-1
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
    
    uopt = [];
    if sigma == 1
        uopt = exp_sol_1_1k(x_hist(:,i)-xs(:,1))+us(1);
    end
    if sigma == 2
        uopt = exp_sol_2_1k(x_hist(:,i)-xs(:,2))+us(2);
    end
    u_hist(:,i) = uopt;
    x_hist(:,i+1) = A*x_hist(:,i) + B*u_hist(:,i);
    
    if(i==33) 
        A = G2.A; B = G2.B; 
    end
    if(i==66)
        A = G1.A; B = G1.B;
    end
    
end
x_hist(:,length(t)) = A*x_hist(:,length(t)-1)+B*u_hist(length(t)-1);



figure;
subplot(3,1,1);
plot(t,x_hist(1,:));
grid on;
ylabel('x_1(mA)');

subplot(3,1,2)
plot(t,x_hist(2,:));
grid on;
ylabel('x_2(V)');

subplot(3,1,3)
plot(t(1:end-1),u_hist);
grid on;
ylabel('Duty cycle(%)');
xlabel('Time(ms)');

figure;
plot(expmpc(1).optimizer.Set);
hold on; grid on;
plot(x_hist(1,:)-xs(1,1),x_hist(2,:)-xs(2,1),'k--o');
ylabel('x_2(V)');
xlabel('x_1(mA)');

