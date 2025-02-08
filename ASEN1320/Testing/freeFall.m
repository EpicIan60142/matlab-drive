clear
clc

h = input("Enter a height in meters:");
time = 0:0.0001:5;
distance = h - 0.5*(9.82)*time.*time;
groundTime = calcHit(distance, time)