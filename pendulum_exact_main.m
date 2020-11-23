clear;close all;clc

N = 30; c = 10; g = 9.81; L = 1; tf = 10;
q11 = 10;
q22 = 10;
th0 = deg2rad(5);      thd0 = 0;
thf = deg2rad(10);               thdf = 0;

t = linspace(0,tf,N);
S = eye(2);
AA = [0   1;
    -g/L 0];
BB = [0; 1];
Q = [q11 0; 0 q22];
R = 1;
tRiccati = tf-t;
[tt,y] = ode45(@(t,y)PendRiccati(y,g,L,AA,BB,Q,R,S),t,[thf; thdf]);



