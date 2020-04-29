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

[K,~,gamma] = mixsyn(ohm100,W1,[],W2);

U = feedback(K,G);
S = feedback(1,G*K);
T = feedback(G*K,1);

figure;
step(U)
figure;
step(T,'r');
K=reduce(K,2);
T = feedback(G*K,1);
hold on;
step(T);
gamma
