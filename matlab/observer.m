clc;
clear;
close all;

% System dynamics
load ohm100_10kHz.mat;
sys = ohm100_10kHz;
A=sys.A;B=sys.B;C=[0,1];
L = place(A',C',[0.5,0.51])'
%L=[0;0.1];
abs(eig(A-L*C))

