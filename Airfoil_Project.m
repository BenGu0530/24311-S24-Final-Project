clc;
clear;
close all;
warning off;

fprintf('\n\n24-311 S24 Final Project');
fprintf('\nBen Gu, Thomas Luo, Henry Perine, Steven Powell\n\n\n');


% Parameters
u0 = 10; % windspeed
u1 = 2*u0/3; % another speed
rho = 10; % some rho


syms R C_L V H; % Length, Coefficient of Lift, Volume, Height
W = V / (R + H);

sum_i = 0;
for i = 1:20
    v = R * i * u0/ 20; % another speed
    A = sqrt(u1^2 + v^2); % first sqrt
    B = 1 / sqrt(1 + (u1 + v)^2); % second sqrt part
    sum_i = sum_i + C_L * rho * A * B * R * W;
end
P = sum_i;

L = @(C_L, R, V) -1 * P; % minimize the negative P
dLdCL = @(C_L, R, V) 0; % 
dLdLength = @(C_L, R, V) 0; % 
dLdVolume = @(C_L, R, V) 0; % 

C_L = 0.5; % Initial Condition in C_L
R = 1; % Initial Condition in R
V = 0.01; % Initial Condition in V

alpha = 0.01;

for iter = 1:1000 
    % Gradient descent
    gradC_L = dLdCL(C_L, R, V);
    gradR= dLdLength(C_L, R, V);
    gradV = dLdVolume(C_L, R, V);
    
    % Update
    C_L = C_L - alpha * gradC_L;
    R = R - alpha * gradR;
    V = V - alpha * gradV;
    
    % Print loss
    currentLoss = L(C_L, R, V);
    fprintf('Iteration %d: Loss = %f\n', iter, currentLoss);
end

fprintf('Optimized CL: %f, Length: %f, Volume: %f\n', C_L, R, V);
