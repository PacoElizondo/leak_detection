%% Nominal Analysis
clear
close all
clc

% Load nominal residuals: complete fault sensitivity matrix (31x31)
load nominal_residuals 

% change these indices when plotting for sensors 10 and 12 instead
sensor1 = 14;
sensor2 = 30;

% Leak sensitivity matrix with selected sensors
Omega = [res_nom(sensor1,:); res_nom(sensor2,:)];

% --- Plot Nominal Residuals ---
figure('Name','Nominal Residuals','Color','w');
scatter(Omega(1,:), Omega(2,:), 50, 'filled', ...
    'MarkerFaceColor', [0.2, 0.45, 0.8]); % use consistent blue tone
hold on
plot(0, 0, 'ko', 'MarkerSize', 8, 'LineWidth', 1.2); % origin
quiver(0, 0, Omega(1,1), Omega(2,1), 0, 'r', 'LineWidth', 1.3); % leak #1
quiver(0, 0, Omega(1,sensor1), Omega(2,sensor1), 0, 'g', 'LineWidth', 1.3); % leak sensor1
quiver(0, 0, Omega(1,sensor2), Omega(2,sensor2), 0, 'm', 'LineWidth', 1.3); % leak sensor2
xlabel('Pressure residual in node 14', 'FontSize', 12);
ylabel('Pressure residual in node 30', 'FontSize', 12);
title('Nominal Residuals for the 31 Leaks', 'FontSize', 14);
legend('Leak residuals','Origin','Leak #1','Leak @sensor1','Leak @sensor2', ...
    'Location','bestoutside');
grid on
box on
axis equal


%% Hanoi Analysis
% f_20 for sensors 14 and 30 
% another f_20 for sensors 10 and 12, the same for f_50
% select desired dataset. Use f_20 for smaller leaks, f_50 for larger ones.

load hanoi_residuals_f_20.mat
% load hanoi_residuals_f_50.mat

% Extract available residuals for the selected sensors
r1 = squeeze(res_dufu(sensor1,:,:));
r2 = squeeze(res_dufu(sensor2,:,:));

% Plot Real Residuals (Hanoi f_20) 
figure('Name','Hanoi Residuals f_20','Color','w');
hold on

N_leaks = 31;
cmap = parula(N_leaks);

% Plot each leak with color and a number label at the trajectory's last point
for leak = 1:N_leaks
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

% Labels & Title
title('Hanoi Residuals $f_{20}$', 'FontSize', 14, 'Interpreter', 'latex');
xlabel(['Pressure residual in node ', num2str(sensor1)], 'FontSize', 12);
ylabel(['Pressure residual in node ', num2str(sensor2)], 'FontSize', 12);

% Minimal legend for just the key reference items
% legend({'Origin', 'Leak #1 direction', ...
%         ['Leak @sensor', num2str(sensor1)], ...
%         ['Leak @sensor', num2str(sensor2)]}, ...
%     'Location', 'northwest');

% Add a colorbar and label leak index
colormap(cmap);
cb = colorbar;
cb.Label.String = 'Leak Number';
cb.Ticks = linspace(0,1,N_leaks);
cb.TickLabels = string(1:N_leaks);

grid on
box on
axis equal

% Zoom: Automatically scales tight to data, with small margins
xrange = [min(r1(:)), max(r1(:))];
yrange = [min(r2(:)), max(r2(:))];
margin_x = 0.05 * range(xrange);
margin_y = 0.05 * range(yrange);
xlim([xrange(1)-margin_x, xrange(2)+margin_x]);
ylim([yrange(1)-margin_y, yrange(2)+margin_y]);



%% Confusion Matrix Computation
% Compute correlation between real residuals and hypotheses
N_residuals = length(r1);
Gamma = zeros(31,31);

for leak = 1:31
    for k = 1:N_residuals
        V_Ro = zeros(31,1);
        for hypothesis = 1:31
            % Normalized correlation between observed and nominal residuals
            V_Ro(hypothesis) = [r1(k,leak), r2(k,leak)] * ...
                [Omega(1,hypothesis), Omega(2,hypothesis)]' / ...
                (norm([r1(k,leak),r2(k,leak)]) * norm([Omega(1,hypothesis),Omega(2,hypothesis)]));
        end
        [~, winner] = max(V_Ro);
        Gamma(leak, winner) = Gamma(leak, winner) + 1;
    end
end


%% Average Topological Distance (ATD)
load matrix_D.mat % D: matrix (31x31) with all node distances
ATD = 0;
atd_vec = zeros(1,31);

for leak = 1:31
    for hypothesis = 1:31
        ATD = ATD + Gamma(leak,hypothesis) * D(leak,hypothesis);
    end
    atd_vec(leak) = ATD / sum(Gamma(leak,:));
end

ATD = ATD / (31 * N_residuals); % According to Remark 3 in report
disp(['Final Average Topological Distance (ATD): ', num2str(ATD,'%.3f')]);
