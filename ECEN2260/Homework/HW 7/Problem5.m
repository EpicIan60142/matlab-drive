clc; clear; close all;

s = tf('s');

w0 = 10;

H1 = 1+(s/w0);
H2 = 1-(s/w0);

figure
hold on
bode(H1)
bode(H2)