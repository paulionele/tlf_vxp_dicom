clear;clc
x = 0:0.01:2*pi;
y = sin(x);

fx = gradient(y,0.01);


plot(x,fx)