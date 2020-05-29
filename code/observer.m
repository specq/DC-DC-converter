clc;
clear;
close all;

% System dynamics
load ohm100_1kHz.mat;
load ohm100_10kHz.mat;
load ohm18_1kHz.mat;
load ohm18_10kHz.mat;
sys = ohm100_10kHz;
A=sys.A;B=sys.B;C=[0,1];
L = place(A',C',[0.5,0.51])'

