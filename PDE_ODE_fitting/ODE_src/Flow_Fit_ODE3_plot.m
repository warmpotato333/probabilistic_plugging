% Aya Yu
% Created Feb 23 2026
% This plot data from Flow_Fit_ODE3_multi

%%%% Process Data %%%%
%% Finding the varience and the mean of the Rmin
 
% Initiate a cell to store all the statistics
ODRstat_all = cell(size(ODRsol_all));

% Iterate through the Rmin_all for Rstat_all
for i = 1:numel(ODRsol_all)
    % Set Rmin to current cell of Rmin_all
    ODRsol = ODRsol_all{i};
    % Find mean
    avg = mean(ODRsol(:, 2, :), 3, 'omitnan');
    % Find population vairence
    varience = var(ODRsol(:, 2, :), 1, 3, 'omitnan');
    % Find standard deviation
    stdDev = std(ODRsol(:, 2, :), 1, 3, 'omitnan');
    % Store statistic into cell array
    ODRstat_all{i} = [avg varience stdDev]; % the order by column is mean, varience, standard deviation
end

% ## Seperate out plugged and unplugged stats ##

% Find indeces with plugged data
plugged_idx = find(ODplugged == 1);
% Find indeces with no plug data
unplugged_idx = find(ODplugged == 0);

% Initiate two cell arrays to store all the statistics (mean, std_dev,
% varience)
ODRstat_plugged = cell(size(ODRstat_all));
ODRstat_unplugged = cell(size(ODRstat_all));

% Iterate through plugged stats
for i = plugged_idx
    ODRstat_plugged{i} = ODRstat_all{i};
end
% Iterate through unplugged stats
for i = unplugged_idx
    ODRstat_unplugged{i} = ODRstat_all{i};
end



%%
%%%% Plot %%%%

%% Plot plugged time series

for i = 1:length(ODRsol_all)
    ODRsol = ODRsol_all{i};
    ODtend = ODtend_all{i};
    plt_time = ODtend;
    % only plot if ODplugged(i) is 1
    if ODplugged(i) == 1
        figure;
        hold on;
        for ii = 1:size(ODRsol, 3)
            %plot(plt_time(1500/h:end), X_vec(1500/h:end, 2, iii));
            plot(plt_time(:), ODRsol(:, 2, ii));
        end
        title('Time Series of plugged data');
        xlabel('Time');
        ylabel('R');
        hold off;
    end
end

%% Plot unplugged time series

for i = 1:length(ODRsol_all)
    ODRsol = ODRsol_all{i};
    ODtend = ODtend_all{i};
    plt_time = ODtend;
    % only plot if ODplugged(i) is 0
    if ODplugged(i) == 0
        figure;
        hold on;
        % iterate through each page of ODRol. This is because state of each
        % wave R_i and x_i are stored in a different page (we are only plotting the R_i in this)
        for ii = 1:size(ODRsol, 3)
            %plot(plt_time(1500/h:end), X_vec(1500/h:end, 2, iii));
            plot(plt_time(:), ODRsol(:, 2, ii));
        end
        title('Time Series of not plugged data');
        xlabel('Time');
        ylabel('R');
        hold off;
    end
end


%% Plot the main, standard deviation for plugged data
for i = 1:numel(ODRstat_plugged)
    ODRstat = ODRstat_plugged{i};
    % If the cell is not empty, plot it.
    if ~isempty(ODRstat) 
        avg = ODRstat(:, 1);
        varience = ODRstat(:, 2);
        stdDev = ODRstat(:,3);
        
        % plot it
        figure;
        hold on;
        plot(avg);
        plot(avg+stdDev);
        plot(avg-stdDev);
        hold off;
    end
end

%% Plot the main, standard deviation for unplugged data
for i = 1:numel(ODRstat_unplugged)
    ODRstat = ODRstat_unplugged{i};
    % If the cell is not empty, plot it.
    if ~isempty(ODRstat) 
        avg = ODRstat(:, 1);
        varience = ODRstat(:, 2);
        stdDev = ODRstat(:,3);
        
        % plot it
        figure;
        hold on;
        plot(avg);
        plot(avg+stdDev);
        plot(avg-stdDev);
        hold off;
    end
end

%% Plot the distribution of the distance function
% figure;
% plt_x = -10:0.01:10;
% plt_y = lambda1 * exp( -lambda2 * (plt_x - 0).^2 );
% plot(plt_x, plt_y);
% title('Distance Function Distribution');

%% Plot phase space
% figure;
% plot3(X_vec(873/h:end, 2, 1), X_vec(873/h:end, 2, 2), X_vec(873/h:end, 2, 3));
% hold on;
% scatter3(X_vec(873/h, 2, 1), X_vec(873/h, 2, 2), X_vec(873/h, 2, 3));   % plot initial condition as a dot
% title('Phase Space, Wave Crest Height VS Speed');
% xlabel('Wave Speed');
% ylabel('Wave Crest Height');
% hold off;
% %% Plot the wave postions
% 
% for i = 1:length(X_all)
%     X_vec = X_all{i};
%     tend = Otend_all{i};
% 
%     figure;
%     subplot(2, 1, 1)
%     hold on;
%     plt_time = tend;
%     for iii = 1:size(X_vec, 3)
%         plot(plt_time, mod(X_vec(:, 1, iii), L));
%     end
%     title('Wave Postion by Time');
%     xlabel('Time');
%     ylabel('Position x');
%     hold off;
% 
%     % Plot the wave crest heights
%     subplot(2, 1, 2);
%     hold on;
%     plt_time = tend;
%     for iii = 1:size(X_vec, 3)
%         %plot(plt_time(873/h:end), x_vec(873/h:end, 2, iii));
%         plot(plt_time(:), X_vec(:, 2, iii));
%     end
%     title('Time Series of R');
%     xlabel('Time');
%     ylabel('R');
%     hold off;
% 
%     toc;
% 
% end