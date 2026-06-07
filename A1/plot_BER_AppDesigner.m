function plot_BER_AppDesigner(app, SUMF, DECO, MMSE, OPT, ref_snr_dB, user_idx)
    % Clear previous plot
    cla(app.UIAxes);

    % Define a color map or use distinguishable colors for each user
    colors = lines(length(user_idx));  % MATLAB's 'lines' colormap to assign different colors
    line_styles = {'-o', '-s', '-^', '-d'};  % Different line styles for SUMF, DECO, MMSE, OPT
    
    % Loop over each user in user_idx and plot the BER for each receiver
    for u_idx = 1:length(user_idx)
        idx = user_idx(u_idx);

        % Check if the user index is valid
        if idx > size(SUMF, 1)
            warning(['User index ', num2str(idx), ' exceeds the number of users. Skipping this user.']);
            continue;
        end

        % Extract BER for the current user
        BER_SUMF = SUMF(idx, :);
        BER_DECO = DECO(idx, :);
        BER_MMSE = MMSE(idx, :);
        BER_OPT  = OPT(idx, :);

        % Plot BER using semilogy in the UIAxes for the current user
        semilogy(app.UIAxes, ref_snr_dB, BER_SUMF, line_styles{1}, 'LineWidth', 1, 'Color', colors(u_idx,:), 'DisplayName', ['User ', num2str(idx), ' - SUMF']);
        hold(app.UIAxes, 'on');
        semilogy(app.UIAxes, ref_snr_dB, BER_DECO, line_styles{2}, 'LineWidth', 1, 'Color', colors(u_idx,:), 'DisplayName', ['User ', num2str(idx), ' - DECO']);
        semilogy(app.UIAxes, ref_snr_dB, BER_MMSE, line_styles{3}, 'LineWidth', 1, 'Color', colors(u_idx,:), 'DisplayName', ['User ', num2str(idx), ' - MMSE']);
        semilogy(app.UIAxes, ref_snr_dB, BER_OPT, line_styles{4}, 'LineWidth', 1, 'Color', colors(u_idx,:), 'DisplayName', ['User ', num2str(idx), ' - OPT']);
    end

    % Add grid and labels in the app UIAxes
    grid(app.UIAxes, 'on');
    xlabel(app.UIAxes, 'SNR (dB)', 'FontSize', 12);
    ylabel(app.UIAxes, 'BER (Bit Error Rate)', 'FontSize', 12);
    title(app.UIAxes, ['BER for Users: ', num2str(user_idx)], 'FontSize', 14);

    % Add legend
    legend(app.UIAxes, 'show', 'Location', 'southwest');
    hold(app.UIAxes, 'off');

    % Set limits for better visualization
    ylim(app.UIAxes, "padded");  % Adjust Y-axis for BER range
    xlim(app.UIAxes, [min(ref_snr_dB) max(ref_snr_dB)]);  % Adjust X-axis to SNR range
    xticks(app.UIAxes, ref_snr_dB)
end
