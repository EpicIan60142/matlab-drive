clc; clear; close all;

s = tf('s');
f0 = 1/(2*pi);
w0 = 2*pi*f0;

A = 10;

% K = (A*1i*w0)/(1+1i);
K = A*w0;

H = (K/s)*(1+(s/w0));

bode(H)

