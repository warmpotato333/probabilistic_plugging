%%% Created on 11/6/2026 %%%
%%% By: Aya Yu %%%
% This code fit the data from PDE simulation with gaussian process
% regression

%Prepare the data
PR_T0 = readtable('~/Documents/GitHub/probabilistic_plugging/GPR/PDE_GPR/PlugRate_Data.csv');
obs_x = [PR_T0.a(446:end)];
obs_y = [PR_T0.plugRate(446:end)];
% fit gpr
gprMdl0 = fitrgp(obs_x, obs_y);

%% compute predicted response and 95% prediction interval
pred_x = linspace(1.31,1.37, 50)';
[pred_y, ~, pred_yint] = predict(gprMdl0, pred_x);

%% plot the prediction
figure(Theme='light');
hold on;
scatter(obs_x, obs_y, 'blue');
plot(pred_x, pred_y, 'black', LineWidth=1);
patch([pred_x;flipud(pred_x)],[pred_yint(:,1);flipud(pred_yint(:,2))],'k','FaceAlpha',0.1);
legend({'Plug rate observations', 'GPR predictions', '95% prediction interval'}, Location='best');
xlabel('a');
ylabel('Plug rate');
xlim([1.31 1.37]);
%hold off;


