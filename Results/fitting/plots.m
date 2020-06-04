clc
clear
close all

load x_hist.mat
load u_hist.mat
load expmpc.mat
load regions.mat

xs = [0.05,5];

x1 = readmatrix('x1.csv');
x2 = readmatrix('x2.csv');
u = readmatrix('u.csv');
reg = readmatrix('regions.csv');

x1(x1(:,1)<0,:) = [];
x2(x2(:,1)<0,:) = [];
u(u(:,1)<0,:) = [];
reg(reg(:,1)<0,:) = [];

x1_d = x1(1,2);
x2_d = x2(1,2);
u_d = [];
reg_d = [];

% Sampling
time = 1000*u(:,1);
T1 = 0.1;
T2 = 0.05;
for i=1:length(time)
    if time(i)>T1
        x1_d = [x1_d;x1(i,2)];
        x2_d = [x2_d;x2(i,2)];
        T1 = T1 + 0.1;
    end
    if time(i)>T2
        u_d = [u_d;u(i,2)];
        reg_d = [reg_d;reg(i,2)];
        T2 = T2 + 0.1;
    end
end
u_d=u_d(1:180,:);

% Computing the regions
reg_d=reg_d(1:180,:);
reg = zeros(180,1);

reg(reg_d<0.25) = 0; 
reg(reg_d>0.25 & reg_d<0.75) = 1;
reg(reg_d>0.75 & reg_d<1.25) = 2;
reg(reg_d>1.25 & reg_d<1.75) = 3;
reg(reg_d>1.75 & reg_d<2.25) = 4;
reg(reg_d>2.25 & reg_d<2.75) = 5;

% Computing x1
x1_d = x1_d(1:180,:); 
x1_d = x1_d/7.5;

x2_d=x2_d(1:180,:);
time = 0:0.1:(length(u_d)-1)*0.1;

figure;
plot(expmpc.optimizer.Set);
hold on;
p1=plot(x1_d-xs(1)*ones(1,length(x1_d)), x2_d-xs(2)*ones(1,length(x2_d)), 'k--o');
p2=plot(x_hist(1,:)-xs(1)*ones(1,180), x_hist(2,:)-xs(2)*ones(1,180), 'b-o');
ylabel('x_2(V)');
xlabel('x_1(mA)');

figure;
subplot(4,1,1);
plot(time,1000*x1_d);
hold on;
plot(time,1000*x_hist(1,:));
grid on;
ylabel('x_1(mA)');


subplot(4,1,2);
plot(time,x2_d);
hold on;
plot(time,x_hist(2,:));
grid on;
ylabel('x_2(V)');

subplot(4,1,3);
plot(time,u_d);
hold on;
plot(time,u_hist);
grid on;
ylabel('Duty cycle(%)');
legend('Real system','Simulation');

subplot(4,1,4);
plot(time,reg);
hold on;
plot(time,regions);
grid on;
ylabel('Region ID');
xlabel('Time(ms)');

    
