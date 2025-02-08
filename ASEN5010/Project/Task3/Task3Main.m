%% ASEN 5010 Project Task 3 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Get [R_sN] and save answer

RsN = calcRsN();

f1 = fopen("RsN_ans.txt", 'w');
ans_RsN = fprintf(f1, "%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f", RsN(1,1), RsN(1,2), RsN(1,3), RsN(2,1), RsN(2,2), RsN(2,3), RsN(3,1), RsN(3,2), RsN(3,3));
fclose(f1);

%% Save answer for omega_{Rs/N}
f2 = fopen("w_ans.txt", 'w');
ans_w = fprintf(f2, "%.7f %.7f %.7f", 0, 0, 0);
fclose(f2)



