%% 1)
clear all
close all
clc
syms x r A
sol=solve(-A*x^2+r*x==0,x)
%% 2)
A=[0 10 50 0 ;...
    -1 3 10 0;...
    -2 10 0 0;...
    0 0 0 1];
r=[5;0;0;1];
x=A\r

