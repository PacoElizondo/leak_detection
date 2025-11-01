clear
close all
clc

load nominal_residuals
load matrix_D.mat

% Define both sensor configurations
configs = {[14, 30], [10, 12]};
config_names = {'Sensors 14 and 30', 'Sensors 10 and 12'};

% Constants

N_NODES = 31;

residual_files = {'hanoi_residuals_f_20.mat', 'hanoi_residuals_f_50.mat'};
residual_names = {'f_20', 'f_50'};
residual_names_title = {'f_{20}', 'f_{50}'};

% Alloc
Gamma = cell(2,2); % Confussion matrices
ATD = zeros(2,2);

for i_config = 1:2
    sensor1 = configs{i_config}(1);
    sensor2 = configs{i_config}(2);

    % Residual vectors matrix nxm [delta_r_i/delta_f_i]_j from i=1:n and
    % j=1:m where n is number of nodes and m is number of leaks (n=m)
    Omega = [res_nom(sensor1,:); res_nom(sensor2,:)]; 

    for i_residuals = 1:2
        load(residual_files{i_residuals})

        r1 = squeeze(res_dufu(sensor1,:,:));
        r2 = squeeze(res_dufu(sensor2,:,:));

        n_residuals = length(r1);
        Gamma_current_config = zeros(N_NODES,N_NODES);

        for leak = 1:N_NODES
            for k = 1:n_residuals
                rho = zeros(N_NODES,1);
                for hypothesis = 1:N_NODES
                    rho(hypothesis) = [r1(k,leak), r2(k,leak)] * ...
                    [Omega(1,hypothesis), Omega(2,hypothesis)]' / ...
                    (norm([r1(k,leak), r2(k,leak)]) * norm([Omega(1,hypothesis), ...
                    Omega(2,hypothesis)]));
                end
                [~, winner] = max(rho);
                Gamma_current_config(leak,winner) = Gamma_current_config(leak,winner) + 1;
            end
        end

        Gamma{i_config,i_residuals} = Gamma_current_config;

        ATD_current_config = 0;
        for leak = 1:N_NODES
            for hypothesis = 1:N_NODES
                ATD_current_config = ATD_current_config + Gamma_current_config(leak, hypothesis) * D(leak,hypothesis);
            end
        end
    
        ATD(i_config, i_residuals) = ATD_current_config / (N_NODES*n_residuals);

    end

end

for i_config = 1:2
    for i_residuals = 1:2
        filename = sprintf('confusion_matrix_sensors_%d_%d_%s.csv', ...
            configs{i_config}(1), configs{i_config}(2), residual_names{i_residuals});
        writematrix(Gamma{i_config, i_residuals}, filename);
    end
end

%% Hardest leaks to separate
% To detect or localize because of sensor location in the network.


fprintf('Most confused leaks (with samples): \n');
for i_config = 1:2
    for data = 1:2
        fprintf('\n%s:\n', config_names{i_config});
        fprintf('%s:\n', residual_names{data});
        
        Gamma_comp = Gamma{i_config, data};

        total_samples = sum(Gamma_comp, 2); % Total samples per leak
        classification_rate = diag(Gamma_comp) ./ total_samples; % Correct localizations over total guesses
        [sorted_rates, sorted_nodes] = sort(classification_rate);

        
        for i = 1:min(10, length(sorted_nodes))
            fprintf('Node %d (%.1f%%) ', sorted_nodes(i), sorted_rates(i)*100);
            fprintf('\n')
        end
        fprintf('\n');
    end
end

%% Plotting


for i_config = 1:2
    sensor1 = configs{i_config}(1);
    sensor2 = configs{i_config}(2);

    Omega = [res_nom(sensor1,:); res_nom(sensor2,:)];

    % Nominal Residuals
    figure('Name','Nominal Residuals','Color','w');
    scatter(Omega(1,:), Omega(2,:), 50, 'filled', ...
        'MarkerFaceColor', [0.2, 0.45, 0.8]); % use consistent blue tone
    hold on
    plot(0, 0, 'ko', 'MarkerSize', 8, 'LineWidth', 1.2); % origin
    quiver(0, 0, Omega(1,1), Omega(2,1), 0, 'r', 'LineWidth', 1.3); % leak #1
    quiver(0, 0, Omega(1,sensor1), Omega(2,sensor1), 0, 'g', 'LineWidth', 1.3); % leak sensor1
    quiver(0, 0, Omega(1,sensor2), Omega(2,sensor2), 0, 'm', 'LineWidth', 1.3); % leak sensor2
    xlabel(sprintf('Pressure residual in node %d', sensor1), 'FontSize', 12);
    ylabel(sprintf('Pressure residual in node %d', sensor2), 'FontSize', 12);
    title('Nominal Residuals for the 31 Leaks', 'FontSize', 14);
    legend('Leak residuals','Origin','Leak at node 1','Leak at sensor 1','Leak at sensor 2', ...
        'Location','bestoutside');
    grid on
    box on
    axis equal

    % Hanoi Residuals
    for i_residuals = 1:2
        load(residual_files{i_residuals})

        r1 = squeeze(res_dufu(sensor1,:,:));
        r2 = squeeze(res_dufu(sensor2,:,:));
        figure();
        hold on

        cmap = parula(N_NODES);

        % Plot each leak with color and a number label at the trajectory's last point
        for leak = 1:N_NODES
            plot(r1(:, leak), r2(:, leak), 'x', ...
                'Color', cmap(leak,:), ...
                'MarkerSize', 6, 'LineWidth', 1);
            % Numbering each leak at last element
            text(r1(end, leak), r2(end, leak), num2str(leak), ...
                'Color', cmap(leak,:), ...
                'FontSize', 8, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        end

        plot(0, 0, 'ko', 'MarkerSize', 8, 'LineWidth', 1.3); % Origin

        % Key arrows for reference leaks
        quiver(0, 0, Omega(1,1), Omega(2,1), 0, 'r', 'LineWidth', 1.3, 'MaxHeadSize', 0.5);       % Leak #1
        quiver(0, 0, Omega(1, sensor1), Omega(2, sensor1), 0, 'g', 'LineWidth', 1.3, 'MaxHeadSize', 0.5); % Leak @sensor1
        quiver(0, 0, Omega(1, sensor2), Omega(2, sensor2), 0, 'm', 'LineWidth', 1.3, 'MaxHeadSize', 0.5); % Leak @sensor2

        title(['Hanoi Residuals' ; string(residual_names_title(i_residuals)) ; string(config_names(i_config))], 'FontSize', 14);
        xlabel(['Pressure residual in node ', num2str(sensor1)], 'FontSize', 12);
        ylabel(['Pressure residual in node ', num2str(sensor2)], 'FontSize', 12);

        margin_x = 0.05 ;
        margin_y = 0.05 ;
        xlim([min(Omega(1, :)) - margin_x, max(Omega(1, :)) + margin_x]);
        ylim([min(Omega(2, :)) - margin_y, max(Omega(2, :)) + margin_y]);

        colormap(cmap);
        cb = colorbar;
        cb.Label.String = 'Leak Number';
        cb.Ticks = linspace(0,1,N_NODES);
        cb.TickLabels = string(1:N_NODES);

        grid on
        box on
        axis equal
    end
end


