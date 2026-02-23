% Aya Yu
% Created Aug 25 2025
% This simulates the ODE system proposed to descrive the PDE flow dynamic in a tube
% This version is based on the previous "_main" but fixed some bugs

%%%% functions %%%%
tic;
% the RHS of ODE
function Rj_1 = f(xj_0, fparams) %Rj_0 is an array of 2 initial conditions, it containes speed and crest height
    % unpack the parameters
    a = fparams{1};
    Rb = fparams{2};
    omega = fparams{3};
    p = fparams{4};
    lambda1 = fparams{5};
    lambda2 = fparams{6};
    Xi = fparams{7};            % xi is a matrix of all other state variables in the form of [xi, Rj; ...]
    dist_func = fparams{8};
    L = fparams{9};             % Tube length
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

    %%%%test%%%%
    A = sqrt(1 - ((1- Rj)^2)/0.25); %sqrt(1 - ((1- Rj)^2)/0.5); % This was changed to create more varience in speed betweeen large and small wave
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
    [scaled_diff, mi]  = min(abs(candidates), [], 3);
    signed_diff = candidates(mi);
    
    if abs(Rj - Ri) > p
        if Rj > Ri 
            sign = 1;
        else
            sign = -1;
        end
    else
        if signed_diff < 0
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
function [x_vec, tend_vec, xp_vec, plugged] = timeSeries(f, inits, h, tmax, plugsens, fparams)     %inits are the initial conditions of state variables in the form of an array [x1 R1; x2 R2; ...]
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
    %Make the second column be random initial conditions for R
    inits(:,2) = (R0_max - R0_min).*rand(initsNum, 1) + R0_min;
end


%% Parameter Settings and Initializations
% ### Since each a values has differnt Rc, this currently only works for
% multi runs with the same a value. ###

% Model Parameters                                                                              % dimensioness film thickness parameter
Bo = 1;                                                                                % Bond number
Rb = 0.8423;                                                                           % bifurcation critical R value
Rc = 0.6691;                                                                           % Critical R which cause the plug
p = 0.001;                                                                             % Dampening term
lambda1 = 0.1;                                                                         % parameter for distance function
lambda2 = 0.5;                                                                         % parameter for disrtance function
la = 12;
L = 2*pi*la;                                                                           % parameter for tube length

% Parameters for a values
numtests=4;                                             % Number of times to run a simulation
atests_min = 1.35;                                      % smallest a value to test
atests_max = 1.35;                                      % largest a value to test
atests_inc = 0.01;                                      % Increment of a value
a_values = atests_min:atests_inc:atests_max;            % create an array of all the a values that will be tested
atests = repelem(a_values, numtests);                   % repeat the a value by numtests times, so this new array could be used by parfor
ODplugged = zeros(1, length(atests));                   % create an array of zeros the same size as atests, when there is a plug, the coursebounding position will be turned to 1, in the end, they are summed up to see how many of each a values are plugged
timeToPlug = nan(1, length(atests));                    % create an array of NaN that marks the the time it take for each test to plug
ODRsol_all = cell(1, length(atests));                   % create a cell that will store all the soluions (Rsol)
ODtend_all = cell(1, length(atests));                   % store tend, the maximal time reached (plugging time) from each run

% Initial Condition Parameters
initsNum = round( L/(2 *sqrt(2)*pi));                            %Number of waves
x0_inc = L/initsNum;                                             %Distance between waves
R0_min = 0.9;                                                    %Min wave R
R0_max = 0.999;                                                  %Max wave R

% Integration Parameters
h = 0.1;                   % dt step
tmax = 4000;                % max integration time
plugsens = 0.7;             % sensitivity of plug detection, how small can R get before simulation stop



%% Run
parfor a_i = 1:length(atests)
    % set a
    a = atests(a_i);                                                                       % a value
    omega = ( (a^2)/(64*Bo * (1-Rc) * (1-Rb) ) ) * (a^4 + 3 - 4*a^2 + 4*log(a));           % Wave growth rate
    fparams = {a, Rb, omega, p, lambda1, lambda2, 0, @distfunc, L, Rc};                    %The position 7 here is set to 0 just as a place holder for xi which would be replaced in @timeSeries()
    % Generate initial conditions
    inits = genInits(initsNum, x0_inc, R0_min, R0_max);
    % Integrate the ODE
    [ODRsol_all{a_i}, ODtend_all{a_i}, xp_vec, plugtrue] = timeSeries(@f, inits, h, tmax, plugsens, fparams);
    % log plugs
    if plugtrue
        ODplugged(a_i) = 1;
    else
        ODplugged(a_i) = 0;
    end
end

% reshape Rsol_all so each row is each individual run, and each column is a different a value
ODRsol_all = reshape(ODRsol_all, numtests, []); 

% process data
plug_rate = sum(ODplugged)/length(ODplugged) % calcuate percentate of runs that resulted in plug


