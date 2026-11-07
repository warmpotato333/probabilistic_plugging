%%% Created on 11/6/2026 %%%
%%% By: Aya Yu %%%
% This code fit the data from PDE simulation with gaussian process
% regression

%Prepare the data
[plugRate, plugConfd] = binofit(numplugs, numtests);   %binofit outputs the plug rate and the confidence interval
% convert to column vector for gpr
obs_x = a_values';    
obs_y = plugRate';
% fit gpr
gprMdl1 = fitrgp(obs_x, obs_y);

%% compute predicted response and 95% prediction interval
pred_x = linspace(1.31,1.37, 50)';
[pred_y, ~, pred_yint] = predict(gprMdl1, pred_x);

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
hold off;


