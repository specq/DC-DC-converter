clc;
close all;
clear;

load ohm18
load ohm100

Ts = ohm18.Ts;

G = stack(1,ohm18,ohm100); 

orderW2 = 6;
[G1u,info1]=ucover(G,ohm18,orderW2);
[G2u,info2]=ucover(G,ohm100,orderW2);
[~,nominal] = min([norm(info1.W1),norm(info2.W1)])

figure;
bodemag(info1.W1,'b',info2.W1,'r');
legend('18 ohms','100 ohms');
W2 = info2.W1;

figure;
bodemag(G/ohm100-1,'b--',W2,'r');
legend('|G/Gn-1|','W2');

s=tf('s');
W1 = c2d((s+10)/2/(s+1e-6),Ts);

% Initialization of the model
P = datadrivenACS;

% Multimodel uncertainty with G1, G2, and G3 
P.Model.Plant = G;

% Frequency vector to be considered
P.Model.Frequency = logspace(-1,log10(pi/Ts),1000);

% Initial controller
order = 2;
epsilon = 1e-3;
z = tf('z',Ts);
Xc = epsilon*z^(order-1);
Yc = z^(order-1)*(z-1);
Kc = Xc/Yc;

% Inialization of the controller
P.setKinit(Kc);
P.Feedback.controller.order = order;

% Objective
P.Feedback.objective.o2W1 = 1/(z-1);

% Constraints
P.Feedback.constraints.cinfW1 = W1;
P.Feedback.constraints.cinfW2 = W2; 
P.Feedback.constraints.cinfW3 = 0;

% Maximum number of iterations
P.Feedback.parameters.maxIter = 100;

% Solve
[Kmat,obj] = solveFB(P);
KFB = computeKFB(P);

% Definition of the feedback transfer functions
Tdd = feedback(G*KFB,1);
Sdd = feedback(1,G*KFB);
Udd = feedback(KFB,G);

opt = stepDataOptions('StepAmplitude',5);
figure;
step(Tdd,opt);

figure;
step(Udd,opt);
