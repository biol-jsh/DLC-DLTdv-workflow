function plot_wing_kinematic(parameter_data, speeds, DLC_color, human_digitizer_color, ylabel_text, current_axes)
    
    axes(current_axes); % Set the current axes
    hold on;
    % Define jitter parameters
    jitter_amount = 0.25;

    % Plot raw data with jitter for automatic treatment
    for i = 1:length(speeds)
        x_jitter = speeds(i) + jitter_amount * (rand(size(parameter_data{i, 1})) - 0.5);
        scatter(x_jitter, parameter_data{i, 1}, 15, 'filled', 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0.5, 'MarkerEdgeColor', DLC_color, 'MarkerFaceColor', DLC_color);
    end

    % Plot raw data with jitter for manual treatment
    for i = 1:length(speeds)
        x_jitter = speeds(i) + jitter_amount * (rand(size(parameter_data{i, 2})) - 0.5);
        scatter(x_jitter, parameter_data{i, 2}, 15, 'filled', 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0.5, 'MarkerEdgeColor', human_digitizer_color, 'MarkerFaceColor', human_digitizer_color);
    end

    % Calculate the mean and 95% confidence interval for the parameter
    mean_parameter = cellfun(@mean, parameter_data);
    mean_mean = mean(mean_parameter(:));
    min_parameter = min(cellfun(@min, parameter_data(:)));
    max_parameter = max(cellfun(@max, parameter_data(:)));
    max_mean_diff = 1.05*max(abs(mean_mean-min_parameter),abs(mean_mean-max_parameter));
    yrange = mean_mean +[-1,1]*max_mean_diff;

    conf_int_parameter = cellfun(@(x) 1.96 * std(x) / sqrt(length(x)), parameter_data);



    % Plot the automatic treatment
    err_auto = errorbar(speeds - 0.02, mean_parameter(:, 1), conf_int_parameter(:, 1), 'DisplayName', 'Automatic', 'Color', DLC_color, 'LineWidth', 1.5);

    % Plot the manual treatment
    err_manual = errorbar(speeds + 0.02, mean_parameter(:, 2), conf_int_parameter(:, 2), 'DisplayName', 'Manual', 'Color', human_digitizer_color, 'LineWidth', 1.5);

    % Customize the plot
    xlabel('Speed (m/s)');
    ylabel(ylabel_text);
    ylim(yrange)
    %title([ylabel_text ' vs. Speed']);
    xlim([2.5, 6.5])
    xticks([3,4.5,6])
    grid on;
    hold on;

    legend([err_auto, err_manual], {'Automatic', 'Manual'});

end
