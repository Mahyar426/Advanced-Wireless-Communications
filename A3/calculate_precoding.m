function [Precoding_Matrix, SINR_dB, P_total_dB, iterations_fixed, iterations_newton, theta, radiation_pattern, elapsed_time] = calculate_precoding(lambda, n_t, H_spacing, alpha, P_noise, d, target_SINR_dB, tolerances)
    tic;
    fixed_iter_tolerance = tolerances(1);
    newton_iter_tolerance = tolerances(2);
    alpha = alpha .* pi ./ 180;
    % Step 1: Build the Channel Matrix H using the given formula
    num_of_users = length(alpha);
    H = zeros(num_of_users, n_t);
    k = 1:n_t;
    amp = (lambda / (sqrt(P_noise) * 4 * pi * d));
    for i = 1:num_of_users
        phase_shift = -1i * 2 * pi * H_spacing * cos(alpha(i));
        H(i, :) = amp .* exp(k .* phase_shift);
    end
    % Step 2: Optimization to Find Precoding Matrix P under SINR Constraints
    gamma_min = 10.^(target_SINR_dB/10); % Convert target SINR to linear scale based on minimum rate requirements
    Precoding_Matrix = zeros(n_t, num_of_users); % Initialize precoding matrix
    iterations_fixed = 0;
    difference = 1;
    lambda_lagrange = ones(2, num_of_users); % Initialize Lagrange multipliers
    % First row is the previous version of lambda as (n)
    % the second row is the one to be updated as (n+1)
    
    while difference > fixed_iter_tolerance % the relative distance by caclulating the norm is measured
        iterations_fixed = iterations_fixed + 1;
        for i = 1:num_of_users
            summation = zeros(n_t);
            for j = 1:length(lambda_lagrange)
                summation = summation + lambda_lagrange(1,j) .* gamma_min .* (H(j,:)' * H(j,:));
            end
            summation = eye(n_t) + summation; % a matrix
            scalar_term = H(i,:) / summation * H(i,:)'; % a number
            lambda_lagrange(2,i) = 1 / ((1 + gamma_min) * scalar_term);
        end
        difference = norm(lambda_lagrange(2,:) - lambda_lagrange(1,:));
        difference = difference / norm(lambda_lagrange(2,:));
        lambda_lagrange(1,:) = lambda_lagrange(2,:);
    end
    % Step3: Moving to Newton
    iterations_newton = 0;
    difference = 1;
    nonlinear_equation = zeros(1, num_of_users);
    nonlinear_equation_diff = zeros(num_of_users);
    while difference > newton_iter_tolerance % all the differences should become below tolerance
        iterations_newton = iterations_newton + 1;
        for i = 1:num_of_users
            summation = zeros(n_t);
            for j = 1:length(lambda_lagrange)
                summation = summation + lambda_lagrange(1,j) .* gamma_min .* (H(j,:)' * H(j,:));
            end
            summation = eye(n_t) + summation; % a matrix
            scalar_term = H(i,:) / summation * H(i,:)'; % a number
            nonlinear_equation(i) = (1 + gamma_min) * scalar_term - 1/lambda_lagrange(1,i);
            nonlinear_equation_diff(i,i) = 1 / (lambda_lagrange(1,i)^2) - (1 + gamma_min) * gamma_min * (scalar_term^2);
            for k = 1: num_of_users
                if i == k
                    continue
                end
                summation = zeros(n_t);
                for j = 1:length(lambda_lagrange)
                    summation = summation + lambda_lagrange(1,j) .* gamma_min .* (H(j,:)' * H(j,:));
                end
                summation = eye(n_t) + summation; % a matrix
                scalar_term = H(i,:) / summation * H(k,:)'; % a number
                nonlinear_equation_diff(i,k) = - (1 + gamma_min) * gamma_min * (scalar_term^2);
            end
            
        end
        update_vector = transpose(nonlinear_equation_diff\transpose(nonlinear_equation)); % inv(d/dx f) * f
        lambda_lagrange(2,:) = lambda_lagrange(1,:) - update_vector;
        difference = norm(update_vector);
        difference = difference / norm(lambda_lagrange(2,:));
        lambda_lagrange(1,:) = lambda_lagrange(2,:);
    end
    % Finally we have lambdas as lambda_lagrange(2,:)
    % Step 4: Calculate W_i and P_i
    summation = zeros(n_t);
    for j = 1:length(lambda_lagrange)
        summation = summation + lambda_lagrange(1,j) .* gamma_min .* (H(j,:)' * H(j,:));
    end
    matrix_in_between = eye(n_t) + summation; % a matrix
    w_coeffs = zeros(n_t, num_of_users);
    A_linear_coeffs = zeros(num_of_users); % matrix of linear equation coefficients
    for i = 1: num_of_users
        for j = 1:num_of_users
            w_coeffs(:,j) = matrix_in_between \ H(j,:)'; % a column vector
            w_coeffs(:,j) = (w_coeffs(:,j))/norm(w_coeffs(:,j)); 
            A_linear_coeffs(i,j) = abs(H(i,:) * w_coeffs(:,j)) ^ 2; % a row vector * column vector = a number
        end
    end
    temp = ones(num_of_users) - eye(num_of_users);
    temp = -temp .* gamma_min;
    temp = temp + eye(num_of_users);
    A_linear_coeffs = A_linear_coeffs .* temp; 
    Powers = linsolve(A_linear_coeffs, ones(num_of_users,1).*gamma_min); % should be a column vectror
    Precoding_Matrix = w_coeffs .* transpose(sqrt(Powers));
    % Up until Powers is 100% correct.
    % Step 5: "Radiation pattern"
    % H_tilda: channel response at each degree
    num_of_points = 2000;
    theta = linspace(0,pi,num_of_points);
    H_tilda = zeros(num_of_points,n_t); 
    gamma_tilda_pol_matrix = zeros(num_of_points,num_of_users);
    phase_shift = -1i * 2 * pi * H_spacing * cos(theta);
    amp = (lambda / (sqrt(P_noise) * 4 * pi * d));
    for k = 1:n_t
        H_tilda(:, k) = amp .* exp(k .* phase_shift);
    end
    gamma_tilda_nominator = abs(H_tilda*Precoding_Matrix) .^ 2; % shape(theta,num_users) matrix
    for i = 1:num_of_users
        gamma_tilda_denominator = sum(gamma_tilda_nominator,2) - gamma_tilda_nominator(:,i) + 1;
        gamma_tilda_pol_matrix(:,i) = gamma_tilda_nominator(:,i) ./ gamma_tilda_denominator; 
    end
    SINR_dB = 10.*log10(max(gamma_tilda_pol_matrix));
    radiation_pattern = gamma_tilda_pol_matrix ./ gamma_min; % For normalization
    %polarplot(theta, radiation_pattern)
    P_total_dB = 10.*log10(sum(Powers));
    elapsed_time = toc;
end