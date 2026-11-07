% Aya Yu
% Created Mar 2 2026
% This plot plug rate and related from data from Flow_Fit_ODE3_multi.m

%% plot the rate of plugging acorss 'a' with error bar
ODnumtest = 100; % only if ODRsol_all is not saved

f = figure;
hold on;
[ODplugRate, ODplugConfd] = binofit(ODnumplugs, ODnumtest);   %binofit outputs the plug rate and the confidence interval

ODplugError = ODplugRate.' - ODplugConfd(:,1);               %these two lines make the confidence interval relative to the center
ODplugError(:, 2) = ODplugConfd(:,2) - ODplugRate.';


h1 = errorbar(ODa_values, ODplugRate, ODplugError(:, 1), ODplugError(:, 2), 'o', 'LineWidth', 0.5, 'Color', 'blue', 'MarkerFaceColor','auto');    %plot with error bar
xlabel('$a$', 'Interpreter', 'latex');
ylabel('Plug formation rate', 'Interpreter','latex');
%title('The Probability of Plug Formation Across different $a$ Values', ['la = ' num2str(la), ', Bo = ', num2str(Bo)], 'Interpreter', 'latex');
f.Theme = 'light';          % light theme for plotting
set(findall(f, '-property', 'FontSize'), 'FontSize', 32);
xlim([min(ODa_values)-0.001, max(ODa_values)+0.001]);
xlim([1.3 1.37]);

% Set the figure name to the title
title_name = ['Plug_Rate_la-', num2str(la), '_Bo-', num2str(Bo)];    
f.Name = title_name;

hold off;

%% Plot histogram of time it took to form a plug

% each cell of ODtend_all don't actually store the last t value, but stores
% the entire time vector. So we stract the last element of each cell, and
% that is the end time where a plug is formed
ODtend_vec = nan(1, length(ODtend_all));
for i = 1:length(ODtend_all)
    % only get the simulations that formed a plug
    if ODplugged(i) == 1 
        ODtend_vec(i) = ODtend_all{i}(end);
    end
end 
% Plot histogram of plug formation times
figure(Theme='light', Position= [100 100 1080 700]);
histogram(ODtend_vec(~isnan(ODtend_vec)), 60);
xlabel('Plug formation time', FontSize=35, FontName = "Times New Romans");
ylabel('Number of sims', FontSize=35, FontName = "Times New Romans");
ax = gca;
ax.FontSize = 30;