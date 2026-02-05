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
Rmin_all = Rmin_all(~cellfun('isempty', Rmin_all));


%% Plot time serieses of wave crests

for col = 1:size(Rmin_all, 2)
    % Get a values
    a = a_values;

    for row = 1:size(Rmin_all, 1)
        figure;
        Rmin = Rmin_all{row, col};                  % Extract Rmin from Rmin_all
        Rtime = 1:size(Rmin, 1);                    % Make time stamps, rn it's just counting from 1, can be changed later
        Rtime = repmat(Rtime.', 1, size(Rmin, 2));  % Repeat the time stamp many times so it can be used to plot in scatter
        
        color = randn(numel(Rmin(:)), 3);          % Random color for each point, *optional*
    
        scatter(Rtime(:), Rmin(:), 5, color, 'filled');

    end
end





%% **For debug, find all the duplicate pairs in Rsol_all**
n = numel(Rsol_all);
pairs = [];

for i = 1:n-1
    for j = i+1:n
        firstRow_i = Rsol_all{i}(1, :);
        firstRow_j = Rsol_all{j}(1, :);

        if isequal(firstRow_i, firstRow_j);   
            [r, c] = ind2sub(size(Rsol_all), [i j]);
            pairs(end+1,:) = [r(1) c(1) r(2) c(2)];
        end
    end
end





