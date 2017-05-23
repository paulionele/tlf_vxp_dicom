clear;clc

myfun = @(x) x.^2;
x = [1,2,3];

x = arrayfun(myfun, x)

x = myfun(x)