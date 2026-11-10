%%% Created on 11/9/2026 %%%
%%% By: Aya Yu %%%
% This code fit the data from PDE simulation with Logistic Regression

%%%% SINGLE REGRESSOR PLUGGED ~ a %%%%
%% Read in data
PR_T0 = readtable('~/Documents/GitHub/probabilistic_plugging/Statistics/PDE_Regression/PDE_Plugged_Rerun1.csv');
PR_T1 =PR_T0(PR_T0.la == 12, :);        % only looking at la= data

%% Fit logistic regression model
T1_model1 = fitglm(PR_T1, "plugged ~ a", Distribution="binomial")

%％ predict and plot
% calculate predictions
pred_x = table(linspace(1.31, 1.37, 50)', VariableNames={'a'});     % X has to be a table
[pred_y, pred_yint] = predict(T1_model1, pred_x, Alpha=0.05);       % prediction

% calculate plug rate of training data
T1_gsum = groupsummary(PR_T1, {'a', 'la'}, 'sum', 'plugged');      % group summary
T1_plugRate = table(T1_gsum.a, T1_gsum.sum_plugged/100, VariableNames={'a', 'plugRate'});   % table for plug rate and a

%% plot
figure(Theme = 'light');
hold on;
plot(pred_x.a, pred_y, Color=[0.8500, 0.3250, 0.0980], LineWidth=1); % plot prediction
patch([pred_x.a(:, 1);flipud(pred_x.a(:, 1))],[pred_yint(:,1);flipud(pred_yint(:,2))],'k','FaceAlpha',0.1, ...
    FaceColor=[0.9290, 0.6940, 0.1250], EdgeColor = [0.9290, 0.6940, 0.1250]); % plot prediction interval
scatter(T1_plugRate.a, T1_plugRate.plugRate);       % plot plug rate from training data
hold off;
xlim([1.31 1.37]);
xlabel('a', 'FontSize', 20);
ylabel('Plug probability', 'FontSize', 20);
legend({ 'Prediction mean','95% prediction interval', 'Plug rate observation'}, Location='best');

%% summary
T1_model1.Rsquared


%%%% MULTI-REGRESSOR MODEL Plugged ~ a + la %%%%