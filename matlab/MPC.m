clc;
clear;
close all;

load ohm100
sys = ohm100;
ref = 5;
A=sys.A; B=sys.B; C = sys.C; D=sys.D;

% Steady state target
ss = [A-eye(2) B;C 0]\[0;0;1]*ref;
xss = ss(1:2);
uss = ss(3);

% From now x and u are relative to xss and uss
% e.g. x <- x-xss and u <- u-uss
% State and input constraints
H = [-1 0;1 0;0 -1;0 1]; h = [0;0.5;0;10]-H*xss;
G = [-1;1]; g=[0;1]-G*uss;

%{
%% Maximal control invariant set

O = polytope(H,h);
figure;
p1 = plot(O,'b');
hold on;

i=1;
while 1
    Oprev = O;
	[F,f] = double(O);	
	% Compute the pre-set
	O = intersect(O, projection(polytope([F*A F*B; zeros(2) G],[f;g]), 1:2));
	if O == Oprev, break; end
	i = i + 1;
end
p2 = plot(O,'r');
legend([p1;p2],{'Feasible set';'Maximal control invriant set'});
%}
%% Terminal set

% Parameters
Q = [10 0; 0 1];
R = 1;
[K,Qf,~] = dlqr(A,B,Q,R);
K=-K;
Ak = A+B*K;

% Computation
figure;
p1 = plot(polytope(H,h),'b');
hold on;
O = polytope([H;G*K],[h;g]);
p2=plot(O,'y');

i=1;
while 1
    Oprev = O;
    [F,f] = double(O);
    O = polytope([F;F*Ak],[f;f]);
    plot(O,'y');
    if O == Oprev
        break;
    end
    i = i+1;
end
p3=plot(O,'r');
legend([p1;p2;p3],{'Feasible set';'Iterations';'Maximal invariant set'});
[Hf,hf] = double(O);

%% Simulation

N = 10;
x = sdpvar(2,N,'full');
u = sdpvar(1,N-1,'full');

obj = u(:,1)'*R*u(:,1);
con = [x(:,2) == A*x(:,1)+B*u(:,1), G*u(:,1)<=g];

for i=2:N-1
    obj = obj + x(:,i)'*Q*x(:,i) + u(:,i)'*R*u(:,i);
    con = [con, x(:,i+1) == A*x(:,i)+B*u(:,i),...
           H*x(:,i)<=h, G*u(:,i) <= g];
end
obj = obj + x(:,N)'*Qf*x(:,N);
con = [con Hf*x(:,N)<=hf];
opt = sdpsettings('solver','gurobi');
ctrl = optimizer(con,obj,opt,x(:,1),u(:,1));

%% Trajectory

x0 = [-0.05;-5];
t = 0:20;
x_hist = zeros(2,length(t)); x_hist(:,1)=x0;
u_hist = zeros(1,length(t)-1);
i = 1;
conv_iter = 0;
while 1
    [uopt,infeasible] = ctrl{x_hist(:,i)};
    if infeasible == 1
        error('Could not solve the problem'); 
    else
        u_hist(i) = uopt;
    end
    x_hist(:,i+1) = A*x_hist(:,i) + B*u_hist(:,i);
    i = i+1;
    if abs(C*x_hist(:,i))<1e-3
        conv_iter = conv_iter + 1;
        if conv_iter == 10
            break;
        end
    else
        conv_iter = 0;
    end
end

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

%{
%% Get some data

Nsamples = 1000;
X = rand(2,Nsamples).*[0.5-xss(1);10-xss(2)]-xss;
U = zeros(1,Nsamples);
for i=1:Nsamples
    [uopt,infeasible]=ctrl{X(:,i)};
    if infeasible == 1
        error('Could not solve the problem');
    else
        U(i) = uopt;
    end
end
%}


