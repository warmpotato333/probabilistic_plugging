% Aya Yu 1/30/2026
% Plotting the data from HPC or local runs

%% Find relative mins - the wave crests

% Create a cell to store all the minimums of R
Rmin_all = cell(size(Rsol_all));

% For loop that start iterate with columns Rsol_all, each column
% coursebounds to a different 'a' value
for col = 1:size(Rsol_all, 2)
    
    % get tas hhs 

    % Get the a value of that column
    a = a_values(col);
    
    % Now iterate through the rows of Rsol_all, each row is a different run
    % with different initial condition but the same 'a' value
    for row = 1:size(Rsol_all, 1)

        % Init several storage variables
        Rsol = Rsol_all{row, col};                      % Extract one Rsol from Rsol_all                     
        Rmat = -1 .* Rsol + a;                          % Turn Rsol into Rmat, which now has bigger R being smaller wave, and smaller R being larger wave
        Rmin = NaN(size(Rmat, 1), size(Rmat, 2));       % Rmin stores the local minima of Rmat, this is to find all the wave crests
        
        % Now iterate each row of Rmat. Each row here represents a
        % snapshot of the fluid at a different time step
        for matrow = 1:size(Rmat, 1)
            
            % Iterate each column within that row. Each column here
            % represents the fluid at different z location
            for matcol = 2:size(Rmat, 2)-1

                % If a point of the fluid is less than its two closests
                % neighbors, save it as a local minimum.
                if Rmat(matrow, matcol) < Rmat(matrow, matcol - 1) && Rmat(matrow, matcol) < Rmat(matrow, matcol + 1)
                    Rmin(matrow, matcol) = Rmat(matrow, matcol);
                end
            end         
        end
        % Store the Rmin matrix into the Rmin_all cell array
        Rmin_all{row, col} = Rmin;
    end
end

% Remove all empty cells (If you only want to look at a portion of the data)
%Rmin_all = Rmin_all(~cellfun('isempty', Rmin_all 3));


%% Plot time serieses of wave crests

for col = 1:size(Rmin_all, 2)
    % Get a values
    a = a_values;

    for row = 1:size(Rmin_all, 1)
        figure;
        Rmin = Rmin_all{row, col};                  % Extract Rmin from Rmin_all
        Rtime = 1:size(Rmin, 1);                    % Make time stamps, rn it's just counting from 1, can be changed later
        Rtime = repmat(Rtime.', 1, size(Rmin, 2));  % Repeat the time stamp many times so it can be used to plot in scatter
        
        color = randn(numel(Rmin(:)), 3);           % Random color for each point, *optional*
    
        scatter(Rtime(:), Rmin(:), 5, color, 'filled');
        %ylim([0.75, 1])

    end
end


%% Find the mean, varience, and std_dev of the Rmin

% Set the filter strength. All local mins with values above this will be
% ignored in stat calculation
filterAbv = 1;

% Initiate a cell to store all the statistics
Rstat_all = cell(size(Rmin_all));

% Iterate through the Rmin_all for Rstat_all
for i = 1:numel(Rmin_all)
    % Set Rmin to current cell of Rmin_all
    Rmin = Rmin_all{i};
    % Filter out values greater then filterAbv by setting them to NaN
    Rmin(Rmin > filterAbv) = NaN;
    % Find mean
    avg = mean(Rmin, 2, 'omitnan');
    % Find population vairence
    varience = var(Rmin, 1, 2, 'omitnan');
    % Find standard deviation
    stdDev = std(Rmin, 1, 2, 'omitnan');
    % Find median
    med = median(Rmin, 2, 'omitnan');
    %Find difference between min and median
    minDiff = abs(med - min(Rmin, [], 2, 'omitnan'));
    % Store statistic into cell array
    Rstat_all{i} = [avg varience stdDev minDiff med]; % the order by column is mean, varience, standard deviation
end

% ## Seperate out plugged and unplugged stats ##

% Find indeces with plugged data
plugged_idx = find(plugged == 1);
% Find indeces with no plug data
unplugged_idx = find(plugged == 0);

% Initiate two cell arrays to store all the statistics (mean, std_dev,
% varience)
Rstat_plugged = cell(size(Rstat_all));
Rstat_unplugged = cell(size(Rstat_all));

% Iterate through plugged stats
for i = plugged_idx
    Rstat_plugged{i} = Rstat_all{i};
end
% Iterate through unplugged stats
for i = unplugged_idx
    Rstat_unplugged{i} = Rstat_all{i};
end

%% Plot the main, standard deviation for plugged data
for i = 1:numel(Rstat_plugged)
    Rstat = Rstat_plugged{i}
    % If the cell is not empty, plot it.
    if ~isempty(Rstat) 
        avg = Rstat(:, 1);
        varience = Rstat(:, 2);
        stdDev = Rstat(:,3);
        
        
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
for i = 201:300 %numel(Rstat_unplugged)
    Rstat = Rstat_unplugged{i}
    % If the cell is not empty, plot it.
    if ~isempty(Rstat) 
        avg = Rstat(:, 1);
        varience = Rstat(:, 2);
        stdDev = Rstat(:,3);
        
        
        % plot it
        figure;
        hold on;
        plot(avg);
        plot(avg+stdDev);
        plot(avg-stdDev);
        hold off;
    end
end







