clc;
clear;
close all;

% System dynamics
load ohm100_10kHz.mat;
sys = ohm100_10kHz;
A=sys.A;B=sys.B;C=[0,1];
L = place(A',C',[0.7,0.8])'

