% Aya Yu
% Created Feb 23 2026
% This plot data from Flow_Fit_ODE3_multi

%%
%%%% Plot %%%%

% % % Plot the distribution of the distance function
figure;
plt_x = -10:0.01:10;
plt_y = lambda1 * exp( -lambda2 * (plt_x - 0).^2 );
plot(plt_x, plt_y);
title('Distance Function Distribution');

%% Plot the wave crest heights

for i = 1:length(ORsol_all)
    X_vec = ORsol_all{i};
    tend = Otend_all{i};
    
    figure;
    hold on;
    plt_time = tend;
    for iii = 1:size(X_vec, 3)
        %plot(plt_time(1500/h:end), X_vec(1500/h:end, 2, iii));
        plot(plt_time(:), X_vec(:, 2, iii));
    end
    title('Time Series of R');
    xlabel('Time');
    ylabel('R');
    hold off;

end

%% Plot phase space
% figure;
% plot3(X_vec(873/h:end, 2, 1), X_vec(873/h:end, 2, 2), X_vec(873/h:end, 2, 3));
% hold on;
% scatter3(X_vec(873/h, 2, 1), X_vec(873/h, 2, 2), X_vec(873/h, 2, 3));   % plot initial condition as a dot
% title('Phase Space, Wave Crest Height VS Speed');
% xlabel('Wave Speed');
% ylabel('Wave Crest Height');
% hold off;
%% Plot the wave postions

for i = 1:length(X_all)
    X_vec = X_all{i};
    tend = Otend_all{i};

    figure;
    subplot(2, 1, 1)
    hold on;
    plt_time = tend;
    for iii = 1:size(X_vec, 3)
        plot(plt_time, mod(X_vec(:, 1, iii), L));
    end
    title('Wave Postion by Time');
    xlabel('Time');
    ylabel('Position x');
    hold off;
    
    % Plot the wave crest heights
    subplot(2, 1, 2);
    hold on;
    plt_time = tend;
    for iii = 1:size(X_vec, 3)
        %plot(plt_time(873/h:end), x_vec(873/h:end, 2, iii));
        plot(plt_time(:), X_vec(:, 2, iii));
    end
    title('Time Series of R');
    xlabel('Time');
    ylabel('R');
    hold off;
    
    toc;

end