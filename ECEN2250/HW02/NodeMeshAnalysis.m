%% ECEN2250 HW02, Ian Faber, 09/06/21
clc; clear; close all;

%% Problem 2: Node Analysis

syms Ib Vx Vs real
syms R1 R2 R3 R4 positive

A = [ (1/R1), ((1/R3) + (1/R4)), (-1/R4);
          -1,                 1,       0;
           0,                 0,       1];
       
b = [Ib;
     Vx;
     Vs];
 
x = A\b;

VA = x(1)
VB = x(2)
VC = x(3)

%% Problem 3: Mesh Analysis

syms Vs Vx real
syms R0 R1 R2 R3 RL positive

A = [ (R0 + R1), -R1, 0;
      -R1, R1 + R2 + R3, -R3;
        0,          -R3, R3 + RL];
    
b = [Vs;
     0;
     -Vx];
 
x = A\b;
 
IA = x(1)
IB = x(2)
IC = x(3)
    

    