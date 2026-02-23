% Aya Yu
% Created Aug 25 2025
% This simulates the ODE system proposed to descrive the PDE flow dynamic in a tube
% This version is based on the previous "_main" but fixed some bugs

%%%% functions %%%%

% the RHS of ODE
function Rj_1 = f(xj_0, fparams) %Rj_0 is an array of 2 initial conditions, it containes speed and crest height
    % unpack the parameters
    a = fparams{1};
    Rb = fparams{2};
    omega = fparams{3};
    p = fparams{4};
    lambda1 = fparams{5};
    lambda2 = fparams{6};
    Xi = fparams{7};            %xi is a matrix of all other state variables in the form of [xi, Rj; ...]
    dist_func = fparams{8};
    L = fparams{9};            % Tube length
    Rc = fparams{10};

    xj = xj_0(1, 1);             
    Rj = xj_0(1, 2);

    Rj_1 = zeros(1,2);           %Initialize an array to store Rj at next time step
    

    % calculate the interaction between Rj and each Ri and sum them up
    total = 0;

    for i = 1:size(Xi,1)
        %product = (Rj - xi(i, 2) + p) * dist_func(xi(i, 1), xj, lambda1, lambda2, la);
        %product = (abs(Rj - Xi(i, 2)) - p)  * dist_func(Xi(i, 1), xj, lambda1, lambda2, L);
        product = abs(abs(Rj - Xi(i, 2)) - p) * sign_func(Xi(i, 1), xj, Xi(i, 2), Rj, p, L) * dist_func(Xi(i, 1), xj, lambda1, lambda2, L);
        total = total + product;
    end
    
    % RHS of ODE, calculate the result at next time step
    A = sqrt(1 - ((1- Rj)^2)/0.5); 
    Rj_1(1,1) = -(1/2) * (1 - (a^2)/A^2 + 2*log(a/A) );
    %Rj_1(1,1) = -(1/2) * (1 - (a^2)/Rj^2 + 2*log(a/Rj) );
    Rj_1(1,2) = ((Rj - 1) * (Rj - Rb) * (Rj - Rc) * omega) - ((Rj -1)*total);  
 

end

% Sign function: determain if a wave should grow or shrink based on their location and size
function sign = sign_func(xi, xj, Ri, Rj, p, L) % Rj is the current wave, Ri is the other wave
    d1 = mod(xj, L) - mod(xi, L);
    d2 = mod(xj, L) - L - mod(xi,L);
    d3 = mod(xj, L) + L - mod(xi, L);
    
    candidates = cat(3, d1, d2, d3);
    abs_diff= min(abs(candidates), [], 3);

    if abs(Rj - Ri) > p
        if Rj > Ri 
            sign = 1;
        else
            sign = -1;
        end
    else
        if xj < xi
            sign = 1;
        else
            sign = -1;
        end
    end
  
end

% Distance function that determains the interaction strength between R_i and R_j
function dist = distfunc(xi, xj, lambda1, lambda2, L)

    d1 = mod(xj, L) - mod(xi, L);
    d2 = mod(xj, L) - L - mod(xi,L);
    d3 = mod(xj, L) + L - mod(xi, L);
    
    candidates = cat(3, d1, d2, d3);
    [scaled_diff, mi]  = min(abs(candidates), [], 3);
    signed_diff = candidates(mi);
    
    dist = lambda1 * exp( -lambda2 * scaled_diff.^2 );
end


% rk4 one step forward
function x_1 = rk4(f, x_0, h, fparams)
    k1 = h * f(x_0, fparams);
    k2 = h * f(x_0 + k1/2, fparams);
    k3 = h * f(x_0 + k2/2, fparams);
    k4 = h * f(x_0 + k3, fparams);
    x_1 = x_0 + (k1 + 2*k2 + 2*k3 + k4)/6;
end

% Integrator, calculate the time series solution from initial conditions
function [x_vec, tend_vec, xp_vec] = timeSeries(f, inits, h, tmax, plugsens, fparams)     %inits are the initial conditions of state variables in the form of an array [x1 R1; x2 R2; ...]
    plugged = false;        % Initialize plugged state
    t_vec = 0:h:tmax;       % Initialize array that store time steps
    ti = 2;                 % initialize time index
    
    inits = permute(inits, [3 2 1]);                     %swap rows with layer
    x_vec = nan(length(t_vec), 2, size(inits, 3));       %initialize tms to store the time series, each layer stores x and R of different waves
    xp_vec = nan(length(t_vec), 2, size(inits,3));       %init vec to store x prime

    x_vec(1, :, :) = inits;                              %set the first row of each layers to the initial condition

    % Calculate the initial rate of change and set that to the first row of
    % xp_vec of each layer

    for i = 1:1:size(inits, 3)
        %put state of all other waves into xi
        xi = squeeze(permute(x_vec(1, :, :), [3 2 1]));         %get all the state at t-1 into xi
        xi(i, :) = [];                                          %take out the state of itself
        fparams{7} = xi;
        
        xp = f(inits(:, :, i), fparams);
        xp_vec(1, :, i) = xp;
    end

    while (ti <= length(t_vec) && plugged == false)
        ii = 1;     %initialize ii counter
        while (ii <= size(x_vec, 3) && plugged == false)
            xj_0 = x_vec(ti-1, :, ii);
            %put state of all other waves into xi
            xi = squeeze(permute(x_vec(ti-1, :, :), [3 2 1]));         %get all the state at t-1 into xi
            xi(ii, :) = [];                                            %take out the state of itself
            fparams{7} = xi;
            
            % update the tms with rk4 
            %disp(xj_0);
            x_vec(ti, :, ii) = rk4(f, xj_0, h, fparams);
            % store x'
            xp_vec(ti, :, ii) = f(xj_0, fparams);
            % Check if the wave is close to plugging, if a plug is about to form, break
            if x_vec(ti, 2, ii) <= plugsens             %This number controls how sensitive is plug detection
                disp(['plugged formed at t = ', num2str(t_vec(ti))]);
                plugged = true;
            end
            ii = ii+1;
        end  
        ti = ti+1;
    end
    % rows with NaN in any column/page for x_vec and xp_vec
    mask_x = any(any(isnan(x_vec),3),2);   
    mask_xp = any(any(isnan(xp_vec),3),2);
    % delete rows with NaN value in the time series for x_vec and xp_vec
    x_vec(mask_x,:,:) = [];               
    xp_vec(mask_xp,:,:) = [];
    tend_vec = t_vec(1:size(x_vec,1));   % return the end ti
end

%Generate Initial Conditions
function inits = genInits(initsNum, x0_inc, R0_min, R0_max)
    inits = zeros(initsNum, 2);     % Make a zero array to stand in for initial conditions

    % Make the first column be wave location spaced by x0_inc for x
    for iii = 2:initsNum
        inits(iii, 1) = inits(iii-1, 1) + x0_inc + 2*randn();
    end
    %M ake the second column be random initial conditions for R
    inits(:,2) = (R0_max - R0_min).*rand(initsNum, 1) + R0_min;
end

%%
%%%% Parameters %%%%

% Model Parameters
a = 1.34;                                                                                % dimensioness film thickness parameter
Bo = 1;                                                                                  % Bond number
Rb = 0.8423;                                                                            % bifurcation critical R value
Rc = 0.6691;                                                                             % Critical R which cause the plug
omega = ( (a^2)/(64*Bo * (1-Rc) * (1-Rb) ) ) * (a^4 + 3 - 4*a^2 + 4*log(a));             % Wave growth rate
p = 0.05;                                                                                 % Dampening term
lambda1 = 0.3;                                                                         % parameter for distance function
lambda2 = 1/pi^2;                                                                          % parameter for disrtance function
la = 12;
L = 2*pi*la;                                                                                  % parameter for tube length



fparams = {a, Rb, omega, p, lambda1, lambda2, 0, @distfunc, L, Rc}; %The position 7 here is set to 0 just as a place holder for xi which would be replaced in @timeSeries()

% Initial Condition Parameters
initsNum = 6; %round( L/(2 *sqrt(2)*pi));                            %Number of waves
x0_inc = L/initsNum;                                             %Distance between waves
R0_min = 0.98;                                                   %Min wave R
R0_max = 1;                                                      %Max wave R

% Integration Parameters
h = 0.1;                   % dt step
tmax = 2000;                % max integration time
plugsens = 0.4;             % sensitivity of plug detection, how small can R get before simulation stop

%%
%%%% Run %%%%

% Generate initial conditions
inits = genInits(initsNum, x0_inc, R0_min, R0_max);
% Integrate the ODE
[X_vec, tend, xp_vec] = timeSeries(@f, inits, h, tmax, plugsens, fparams);




%%
%%%% Plot %%%%

% % Plot the distribution of the distance function
% figure;
% plt_x = -10:0.01:10;
% plt_y = lambda1 * exp( -lambda2 * (plt_x - 0).^2 );
% plot(plt_x, plt_y);
% title('Distance Function Distribution');

%% Plot the wave crest heights
figure;
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
