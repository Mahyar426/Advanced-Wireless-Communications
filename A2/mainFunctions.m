function [Capacity_Wo, Capacity_W, Outage_Capacity_W_CSIT, Outage_Capacity_Wo_CSIT, Egodic_Capacity_W_CSIT, Egodic_Capacity_Wo_CSIT] = mainFunctions(n_t, n_r, SNR_dB, seed, N_samples, perc_cap)   
    % Convert dB to linear scale
    SNR = 10^(SNR_dB / 10);


    rng(seed);

    % Channel matrices H 
    H = zeros(n_r, n_t, N_samples);
    for i = 1:N_samples
        H(:, :, i) = (randn(n_r, n_t) + 1i * randn(n_r, n_t)) / sqrt(2);
    end

    Capacity_Wo = zeros(1, N_samples);
    Capacity_W = zeros(1, N_samples);

    for i = 1:N_samples
        H_herm = H(:, :, i)';
        HH_herm = H(:, :, i) * H_herm;

        % Capacity without CSIT
        Capacity_Wo(i) = real(log2(det(eye(n_r) + (SNR / n_t) * HH_herm)));

        % Capacity with CSIT 
         [~, S, ~] = svd(H(:, :, i));
         lambda_vals = diag(S).^2;
         gamma = water_filling(lambda_vals, SNR);
         capacity_CSIT = sum(log2(1 + gamma .* lambda_vals));
         Capacity_W(i) = capacity_CSIT;



    end

    Capacity_Wo = sort(Capacity_Wo);
    Capacity_W = sort(Capacity_W);

    %  Empirical cumulative distribution functions
    ecdf_values = (1:N_samples) / N_samples;

    % Outage capacity
    Outage_Capacity_Wo_CSIT = interp1(ecdf_values, Capacity_Wo, perc_cap / 100);
    Outage_Capacity_W_CSIT = interp1(ecdf_values, Capacity_W, perc_cap / 100);

    %Ergodic capacity
    Egodic_Capacity_Wo_CSIT = mean(Capacity_Wo);
    Egodic_Capacity_W_CSIT = mean(Capacity_W);

end


%% Water-Filling Function
function gamma = water_filling(lambda_vals, P_t)
    n = length(lambda_vals);
    lambda_vals = sort(lambda_vals, 'descend');
    gamma = zeros(n,1);
    epsilon = 1e-6;
    max_iterations = 1000;
    mu_lower = 0;
    mu_upper = P_t + 1 / min(lambda_vals);
    
    for iter = 1:max_iterations
        mu = (mu_lower + mu_upper) / 2;
        gamma = max(mu - 1 ./ lambda_vals, 0);
        total_power = sum(gamma);
        if abs(total_power - P_t) < epsilon
            break;
        elseif total_power > P_t
            mu_upper = mu;
        else
            mu_lower = mu;
        end
    end
end

