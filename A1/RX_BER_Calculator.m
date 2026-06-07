function [SUMF, DECO, MMSE, OPT, Near_Far] = RX_BER_Calculator(Seq_seed, Noise_seed, seq_len, A_k_db, symbol_num, ref_snr_dB)
    profile on
    sim_iteration = 1;
    A_k = 10.^(A_k_db./20);
    A_k_matrix = diag(A_k);
    ref_snr = 10.^(ref_snr_dB./10);
    sigma = 1./sqrt(ref_snr);
    num_users = length(A_k);

    SUMF = zeros(num_users,length(ref_snr_dB));
    DECO = zeros(num_users,length(ref_snr_dB));
    MMSE = zeros(num_users,length(ref_snr_dB));
    OPT  = zeros(num_users,length(ref_snr_dB));
    
    % Seq generation
    rng(Seq_seed);
    sig_wave = (-1).^randi(2,num_users,seq_len);

    % Correlation matrix
    R = sig_wave*sig_wave'/seq_len;
    R_inv = inv(R);
    [Q,lambda] = eig(R);
    R_square_root = Q * (lambda^(1/2)) * transpose(Q);
    R_square_root_inv = R_square_root^-1;

    % Near Far Resistance
    Near_Far = 1./diag(R_inv);

    % possibilities
    b_k_all_pos = 1-2*de2bi(0:2^num_users-1)';

    % BER simulation
    rng(Noise_seed); % Noise seed
    for i = 1:length(ref_snr)
         % for j = 1:sim_iteration

        % Bit sequence
        b_k =(-1).^randi(2,[num_users,symbol_num]);
        
        % White Gaussian Noise Generation
        Z_tilda = randn(num_users, symbol_num);

        % Correlated noise generation
        sigma_correlated = sigma(i).*R_square_root;
        Z = sigma_correlated * Z_tilda;

        % RX input
        y = R * (A_k_matrix * b_k) + Z;

        % SUMF
        Received_Seq = sign(y);
        SUMF(:,i) = SUMF(:,i) + sum(Received_Seq~=b_k,2);

        % OPT
        Received_Seq = zeros(num_users,symbol_num);
        for k=1:symbol_num
            distance = sum((R_square_root_inv * y(:,k) - R_square_root * A_k_matrix * b_k_all_pos).^2,1);
            [~,ML_arg] = min(distance);
            Received_Seq(:,k) = b_k_all_pos(:, ML_arg);
        end
        OPT(:,i) = OPT(:,i) + sum(Received_Seq~=b_k,2);

        % DECO
        Received_Seq = sign(R_inv * y); 
        DECO(:,i) = DECO(:,i) + sum(Received_Seq~=b_k,2);

        % MMSE
        F = (A_k_matrix^-1) * ((R+(sigma(i)^2)*(A_k_matrix^-2))^-1);
        Received_Seq = sign(F*y);
        MMSE(:,i) = MMSE(:,i) + sum(Received_Seq~=b_k,2);

         % end
    end
        
  
    % BER Calculation 
    SUMF = SUMF./ (symbol_num*sim_iteration);
    DECO = DECO./ (symbol_num*sim_iteration);
    MMSE = MMSE./ (symbol_num*sim_iteration);
    OPT  = OPT./ (symbol_num*sim_iteration);
    profile off

end