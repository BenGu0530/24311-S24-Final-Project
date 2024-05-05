clc;
close all;
clear;


initial_guess = [0.3; 10; 0.084]; % Initial guess for C_L, R, H
max_iters = 100000;
lr_CL = 0.0001; % Learning rate for C_L
lr_R = 0.01;  % Learning rate for R
lr_H = 0.000000001;  % Learning rate for H
max_volume = 100; % Example maximum volume constraint

[x_opt, max_power] = gradientAscentOptimizer(initial_guess, max_iters, lr_CL, lr_R, lr_H, max_volume);
fprintf('Maximum power achieved: %f with parameters [C_L, R, H, V] = [%f, %f, %f, %f]\n', max_power, x_opt(1), x_opt(2), x_opt(3), x_opt(4));

