% Aya Yu 1/30/2026
% Plotting the data from HPC or local runs
% Focuses on analyzing plug rate


%% plot the rate of plugging acorss 'a' with error bar
f = figure;
hold on;
[plugRate, plugConfd] = binofit(numplugs, numtests);   %binofit outputs the plug rate and the confidence interval
plugError = plugRate.' - plugConfd(:,1);               %these two lines make the confidence interval relative to the center
plugError(:, 2) = plugConfd(:,2) - plugRate.';
h1 = errorbar(a_values, plugRate, plugError(:, 1), plugError(:, 2), 'o', 'LineWidth', 0.5, 'Color', '#C7AD16');    %plot with error bar
xlabel('$a$', 'Interpreter', 'latex');
ylabel('Plug formation rate', 'Interpreter','latex');
%title('The Probability of Plug Formation Across different $a$ Values', ['la = ' num2str(la), ', Bo = ', num2str(Bo)], 'Interpreter', 'latex');
f.Theme = 'light';          % light theme for plotting
set(findall(f, '-property', 'FontSize'), 'FontSize', 32);
xlim([min(a_values), max(a_values)]);
%xlim([1.14 1.5]);

% Set the figure name to the title
title_name = ['Plug_Rate_la-', num2str(la), '_Bo-', num2str(Bo)];    
f.Name = title_name;

% %%
% legend([h1, h2, h3, h4], {'Bo = 0.25', 'Bo = 0.5', 'Bo = 1', 'Bo = 2'}, 'FontSize', 30);

%% plot histogram of tipping tipping for a given slice of 'a'
% ## NEED FIX ##
figure;
hist_aind = 2;         %index of 'a' slice
timeToPlug_re = reshape(timeToPlug, numtests, []).';    %reshape timeToPlug so each row is the time to plug for a give 'a'
histogram(timeToPlug_re(hist_aind, :));                    %plot histogram
title('frequency of time to tipping for', ['a = ', num2str(a_values(hist_aind)), ', la = ', num2str(la)]);
xlabel('Tippnig Time');
ylabel('Frequency');


%% Plot histogram of tipping time for a range of 'a'
figure('Theme', 'light');
hist_amin = 15;
hist_amax = 70;
timeToPlug_re = reshape(timeToPlug, numtests, []).'; %reshape timeToPlug so each row is the time to plug for a give 'a'
histogram(timeToPlug_re(hist_amin:hist_amax, :), 60, 'FaceColor', 'k');
%title('frequency of time to tipping for', ['a = ', num2str(a_values(hist_amin)), ' to ', num2str(a_values(hist_amax)), ', la = ', num2str(la)]);
xlabel('Plug formation time', 'FontSize', 20, 'FontName', 'Times New Roman');
ylabel('Frequency', 'FontSize', 20, 'FontName', 'Times New Roman');

%% Plot mean and standard deviation for tipping time at each 'a'

% Compute Statistics
timeToPlug_re = reshape(timeToPlug, numtests, []).';     %reshape timeToPlug so each row is the time to plug for a give 'a'
meanTipTimes = nanmean(timeToPlug_re, 2);
medianTipTimes = nanmedian(timeToPlug_re, 2);
stdTipTimes = nanstd(timeToPlug_re, 0, 2);

%plot
f = figure;
hold on;
% Mean with error bars (std deviation)
errorbar(a_values, meanTipTimes, stdTipTimes, 'o', 'DisplayName', 'Mean ± Std');
%Median
plot(a_values, medianTipTimes, 's', 'LineWidth', 1.5, 'DisplayName', 'Median')

legend show;
xlabel('a');
ylabel('Tipping Time');
title('Tipping Time Across a');
hold off;

% Set the figure name to the title
title_name = gca().Title.String;    
f.Name = title_name;

%% plot box and wisker plot for tipping time at each 'a'
f = figure;

% Plot Box Plot
b = boxchart(atests, timeToPlug);

% Set box color to purple
b.BoxFaceColor = [0.5 0 0.5];  % RGB for purple

% Customize outliers
b.MarkerStyle = 'x';           % x shape
b.MarkerColor = [1 0.5 0];     % orange (RGB)

xlabel('a');
ylabel('Tipping Time');
title('Tipping Time Across a');

% Set the figure name to the title
title_name = gca().Title.String;    
f.Name = title_name;
