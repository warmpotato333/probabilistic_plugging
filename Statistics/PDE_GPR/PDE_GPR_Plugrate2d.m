%%% Created on 11/6/2026 %%%
%%% By: Aya Yu %%%
% This code fit the data from PDE simulation with gaussian process
% regression
% It takes two regressors

% Prepare the data
PR_T0 = readtable('~/Documents/GitHub/probabilistic_plugging/Statistics/PDE_GPR/PlugRate_Data.csv');
PR_T1 = PR_T0(1:5:end, :);  % remove some rows to prevent overfitting
obs_x = [PR_T1.a, PR_T1.la];
obs_y = [PR_T1.plugRate];


%% Fit gpr
gprMdl2 = fitrgp(obs_x, obs_y, KernelFunction='ardsquaredexponential', ...
    Standardize=true, ...
    KernelParameters = []);


%% Affix la, compute and plot prediction along a, 
plt_la = 12;
pred_x1 = [linspace(1.31, 1.37, 50)', repelem(plt_la, 50)'];
[pred_y1, ~, pred_y1int] = predict(gprMdl2, pred_x1);
figure(Theme='light');
hold on
% plot
plot(pred_x1(:, 1), pred_y1, '-', LineWidth=1, ...
    Color=[0.8500, 0.3250, 0.0980]);                                                         %plot prediction
scatter( PR_T1.a(PR_T1.la == plt_la), ( PR_T1.plugRate( PR_T1.la == plt_la ) ), ...
    MarkerEdgeColor=[0.6350, 0.0780, 0.1840]);                                               % plot training observations
patch([pred_x1(:, 1);flipud(pred_x1(:, 1))],[pred_y1int(:,1);flipud(pred_y1int(:,2))],'k','FaceAlpha',0.1, ...
    FaceColor=[0.9290, 0.6940, 0.1250], EdgeColor = [0.9290, 0.6940, 0.1250]);               % plot 95% prediction interval
% styling
legend({ 'GPR predictions', 'Plug rate observations','95% prediction interval'}, Location='best');
xlabel("a", FontSize=20);
ylabel('Plug rate', FontSize=20);
xlim([1.31 1.37]);
t = title(['la = ', num2str(plt_la)], fontsize = 20);
t.HorizontalAlignment = 'left'
t.Units = 'normalized';
t.Position(1) = 0;
hold off

%% Affix a, compute and plot prediction along la
plt_a = 1.35;
pred_x2 = [repelem(plt_a, 50)', linspace(6,12,50)'];
[pred_y2, ~, pred_y2int] = predict(gprMdl2, pred_x2);
figure(Theme='light');
plot(pred_x2(:, 2), pred_y2, 'k');
xlabel('la', FontSize = 20)
ylabel('Plug rate', FontSize = 20)
t = title(['a = ', num2str(plt_a)], fontsize = 20);
t.HorizontalAlignment = 'left'
t.Units = 'normalized';
t.Position(1) = 0;

%% compute and plot prediction surface
[LA, A] = meshgrid(linspace(1.31, 1.37, 50), linspace(4, 12, 50));
pred_x = [LA(:), A(:)];
[pred_y, ~, pred_yint] = predict(gprMdl2, pred_x);
PRED_Y = reshape(pred_y, size(LA));
figure(Theme='light');
hold on;
scatter3(PR_T0.a, PR_T0.la, PR_T0.plugRate, 5, 'blue');
surf(LA, A, PRED_Y);
view(3);
grid on;
hold off;









%%
% plot(pred_x, pred_y, 'black', LineWidth=1);
% patch([pred_x;flipud(pred_x)],[pred_yint(:,1);flipud(pred_yint(:,2))],'k','FaceAlpha',0.1);
% legend({'Plug rate observations', 'GPR predictions', '95% prediction interval'}, Location='best');
% xlabel('a');
% ylabel('Plug rate');
% xlim([1.31 1.37]);
% hold off;
% 
