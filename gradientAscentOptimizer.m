function [x_opt, fval] = gradientAscentOptimizer(initial_guess, max_iters, lr_CL, lr_R, lr_H, max_volume)
    % Parameters
    u0 = 10; % Wind speed
    u1 = 2 * u0 / 3; % Maximum speed due to Betz limit
    rho = 1.204; % Air density (kg/m^3)
    alpha = 7.4; % Reduction factor for top length of the trapezoid per unit height regarding to the height increment

    % Unpack initial_guess
    C_L = initial_guess(1);
    R = initial_guess(2);
    H = initial_guess(3);
    
    %Initialize the basic geometry constraints
    W = 2.8; % standard width (meter)
    B = 1.05 * W; % Cross section bottom width (meter)
    T = 0.5 * W; % Cross section top width (meter)

    H_min = 0.03 * W; % Cross section height min (meter)
    H_max = 0.20 * W; % Cross section height max (meter)

    R_max = 150; % Length max (meter)
    R_min = 0; % Length min (meter)

    CL_max = 1.6; % Coefficient of lift max
    CL_min = 0.2; % Coefficient of lift min
    
    max_power = 0;
    best_params = initial_guess;

    % for plotting
    plt_iter = [];
    plt_power = [];


    for iter = 1:max_iters
        % Calculate power and gradient
        [power_val, grad] = powerFunction(C_L, R, H, u0, u1, rho);

        % Initial parameter update in the direction of the gradient
        C_L_new = C_L + lr_CL * grad(1);
        R_new = R + lr_R * grad(2);
        H_new = H + lr_H * grad(3);

        % Compute top length T and volume V
        dH = abs(H_new - H);
        T = max(0, T - alpha * dH);
        V = (B + T) * H_new / 2 * R_new; % Volume of the trapezoid

        % Apply constraints
        if V > max_volume
            gradRatio = abs(grad(2)) / abs(grad(3));  % Ratio of gradients for R and H

            if gradRatio > 1
                % R gradient dominates; decrease H more to allow increase in R
                adjustmentFactor = lr_R * (1 + gradRatio);  % Scale adjustment based on dominance
                H  = max(H_min, H - lr_H * adjustmentFactor * grad(3));  % More reduction in H
                R = min(R_max, R + lr_R * lr_R * grad(2));  % Potential increase in R
            else
                % H gradient dominates; decrease R more to allow increase in H
                adjustmentFactor = lr_H * (1 + 1 / gradRatio);  % Scale adjustment based on dominance
                R = max(R_min, R - lr_R * adjustmentFactor * grad(2));  % More reduction in R
                H = min(H_max, H + lr_H * lr_H * grad(3));  % Potential increase in H
            end
        % If volume is within limits, update parameters
        elseif V <= max_volume
            plt_iter = [plt_iter;iter];
            plt_power = [plt_power;power_val];
            if max_volume - V > 0 && max_volume - V < 0.05
                gradRatio = abs(grad(2)) / abs(grad(3));
                if gradRatio == 1
                    max_power = power_val;
                    best_params = [C_L; R; H; V];
                    x_opt = best_params;
                    fval = max_power; 
                    return;
                end
            end

            C_L = max(CL_min, min(CL_max, C_L_new));  
            R = max(R_min, min(R_max, R_new));        
            H = max(H_min, min(H_max, H_new));    
            
            if power_val > max_power
                max_power = power_val;
                best_params = [C_L; R; H; V];
            end

        end

        % Display current iteration info
        fprintf('Iteration %d: Power = %f, C_L = %f, R = %f, H = %f, V = %f, T = %f\n', ...
                iter, power_val, C_L, R, H, V, T);
    end

    % Optimal values
    x_opt = best_params;
    fval = max_power; % Final value of the objective function

    % Plot power
    figure()
    plot(plt_iter,plt_power,'b')
    hold on
    yline(max_power,'r')
    hold off
    xlim([0,3000]);
    xlabel('Iteration')
    ylabel('Power (W)')
    title('Power of one blade over iterations')
    legend('Power Over Iterations','Max Power')
end

function [power_val, grad] = powerFunction(C_L, R, H, u0, u1, rho)
    sum_i = 0;
    grad_sum = zeros(3, 1);  % Gradient for C_L, R, H
    for i = 1:20
        v = R * i * u0 / 20;
        A = sqrt(u1^2 + v^2);
        B = 1 / sqrt(1 + (u1 / v)^2);
        sum_component = C_L * rho * A * B * R * H;
        sum_i = sum_i + sum_component;

        % Gradient calculation for each component
        dPdCL = rho * A * B * R * H;
        dPdR = C_L * rho * (A * B * H + A * B * R * (i * u0 / 20) / A * H - sum_component / R);
        dPdH = C_L * rho * A * B * R ;

        % Accumulate gradients
        grad_sum(1) = grad_sum(1) + dPdCL;
        grad_sum(2) = grad_sum(2) + dPdR;
        grad_sum(3) = grad_sum(3) + dPdH;
    end

    power_val = sum_i;
    grad = grad_sum; % Return the gradient vector
end
