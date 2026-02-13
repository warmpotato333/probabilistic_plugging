% Reed Ogrosky, Aya Yu
% July 16th 2025 

tic;
clear ProgressUpdate;

% user = getenv('USER');
% node = getenv('HOSTNAME');
% parpool("Processes");

numtests=10;                                           % Number of times to run a simulation
atests_min = 1.35;                                      %smallest a value to test
atests_max = 1.36;                                      %largest a value to test
atests_inc = 0.005;                                    %Increment of a value
a_values = atests_min:atests_inc:atests_max;            %create an array of all the a values that will be tested
atests = repelem(a_values, numtests);                   %repeat the a value by numtests times, so this new array could be used by parfor
plugged = zeros(1, length(atests));                     %create an array of zeros the same size as atests, when there is a plug, the coursebounding position will be turned to 1, in the end, they are summed up to see how many of each a values are plugged
timeToPlug = nan(1, length(atests));                    %create an array of NaN that marks the the time it take for each test to plug
Rsol_all = cell(1, length(atests));                     %create a cell that will store all the soluions (Rsol)

dq = parallel.pool.DataQueue;                           %start a data queue to keep progress
afterEach(dq, @(~) ProgressUpdate(length(atests)));     %When dq is sent, call ProgressUpdate()


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_plots =400;             % Number of solution snapshots to save
la = 12;
Bo=1;
numptsper2pi = 32;          % Set spatial resolution - number of z gridpoints per 2pi
Ltilde=0;                   %Slip Length for Slip BC Model this is Lambda tilde

% Generate a random seed for the randomStream that would be used in parfor
rng("shuffle");
randSeed = randi(1e9); % randSeed saved for debug

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parfor ww=1:length(atests)

    % generate the random stream, and assign a substream for this worker
    % this should prevent any possilble duplication
    randomStream = RandStream("mrg32k3a", "Seed", randSeed);
    randomStream.Substream = ww; 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a=atests(ww);   % Film thickness ratio (tube radius to core radius) - must be real number greater than 1
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Set spatial resolution
    L = 2*pi*la;              % Don't change - Domain is la*(2*pi)
    numptsper2pi = 32;        % Set spatial resolution - number of z gridpoints per 2pi
    N = numptsper2pi*la;      % Don't change - Total number of z gridpoints over entire domain
    dz = L/N;                 % Don't change - Space between successive z gridpoints
    z = 0:dz:L-dz;            % Don't change - Vector containing all z gridpoints
    num_refinements=0;
    
    Rsolmult=zeros(num_plots+1,N,numtests); 
    wave_crests_mult=zeros(num_plots+1,round(la/sqrt(2)),numtests);
    wave_crest_locs_mult=zeros(num_plots+1,round(la/sqrt(2)),numtests);
    mean_thicknesses_mult=zeros(num_plots+1,round(la/sqrt(2)),numtests);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    disp('***************')
    disp(strcat(['Test #',int2str(ww),'/',int2str(numtests)]));
    disp('***************')


    % Set initial temporal resolution 
    t0 = 0;                                  % Don't change - Start time
    tend=10; %250*round(1/(a-1)^3);       % End time
    tnumpersec1=ceil(8000*(a-1)^3);          % Number of time steps per time unit (initially)
    tnum = round(tnumpersec1*(tend-t0));     % Don't change - Number of time intervals
    dt = (tend-t0)/tnum;                     % Don't change - Time-step
    time_snaps=[0:tend/num_plots:tend];      % Don't change - Vector of times when solution is saved
    t=t0;                                    % Don't change - Set t=initial time
    num_time_steps = 0;                      % Don't change - Keep track of how many time steps taken
    num_dt_halves = 0;                       % Don't change - How many times have we had to halve dt? 
    num_dt_halves_max = 8;                   % How many times are we willing to halve dt?
    
    tvec=[];
    tvec(1)=t0;
    dtsave=tend/num_plots;
    dthalvinglastr=dt;
    
    
    % Don't need to change this
    rand_num_vec=[0.6470 0.7988 0.4049 0.4824 0.0639... % vector of 'random' numbers
        0.8029 0.3511 0.8493 0.9918 0.4843 0.3564...
        0.0544 0.9141 0.7050 0.4381 0.9681 0.7238...
        0.6568 0.5148 0.8771 0.0602 0.4715 0.3110...
        0.1947 0.1510 0.9014 0.0707 0.9584 0.6999...
        0.0120 0.0049 0.8533 0.4673 0.8367 0.2206...
        0.8060 0.2486 0.9983 0.9379 0.5540 0.9127...
        0.6355 0.7801 0.5672 0.1918 0.9062 0.5149];
    
    %Construct initial condition - R=1+some small perturbations
    beta=0.01;               % amplitude of small perturbations
    num_modes=5*la;             % Number of Fourier modes to add to initial condition
    R=1;                     % Don't change - Mean value of R
    for m=1:1:num_modes
        R=R+(a-1)*beta*rand*cos(m*z/la+10000*rand*pi/2);
    end
    
    %Initial condition to check linear stability results
    %R=1+0.01*cos(z/2)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% Probably can keep everything here fixed for now %%%%%%%%%%
    % Create filename to use for saving results - don't need to change--Mark
    % changed to include t
    fname = strcat(['Rsol_Bo',num2str(Bo),'_a',num2str(a),'_la',num2str(la),'_tend',num2str(tend)]);
    fname=strrep(fname,'.','p');   % Replace decimals with 'p' (for point)
    fname=strrep(fname,'-','n');   % Replace negative signs with 'n'
    
    % Don't need to change these
    da_fact=2.05;                           % De-aliasing factor - keep at something larger than 2
    num_times_after_wave_tracking = 50;     % number of times after beginning tracking
    num_times_before_wave_tracking = 4;     % num_times_before_wave_tracking_mat(stcol);  %number of times before beginning tracking
    coarse_factor = 2;                      % Keep fixed at 2 - when necessary, cuts dt by factor of 2
    r = 1;                                  % Which snapshot did we last save?
    
    % Initialize solution variables - Don't need to change these
    Rsol = zeros(num_plots+1,N);   % Solution is saved in Rsol at num_plots times
    Rsol(1,:) = a-R(:);            % First saved solution is just initial condition
    
    % Initialize flux, mean_R variables - keep these fixed for now
    mean_R = zeros(1,num_plots+1);                  % Saves mean value of R at each snapshot
    liq_volume_temp = sum(pi*(a^2-R.^2)*dz);        % Compute film volume over domain
    mean_R_temp_sq= (-liq_volume_temp/(pi*L)+a^2);  % Compute mean of R^2
    mean_R(1,1) = sqrt(mean_R_temp_sq);             % Compute mean of R
    
    % Initialize Noise-Induced Tipping variables
    wave_crests = zeros(num_plots+1,round(la/sqrt(2)));
    wave_crest_locs = zeros(num_plots+1,round(la/sqrt(2)));
    mean_thicknesses = zeros(num_plots+1,round(la/sqrt(2)));
    
    % Initialize wave speed variables - keep these fixed for now
    total_wave_dist = 0; 
    total_distance_traveled = 0;
    wave_speed_start_r = -1;
    wave_speed_start_time = 0;
    [m11 temp1_max] = min(R(:));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Print some initial information to the screen
    disp('-----------------------------------------------');
    disp(strcat(['Solving eq. (13) from Camassa et al. (2012)']));
    disp(strcat(['a=',num2str(a),', Bo=',num2str(Bo),', la=',num2str(la)]));
    disp(strcat(['dt=',num2str(dt),', End time=',num2str(tend),', r_max=',num2str(num_plots)]));
    disp('-----------------------------------------------');
    disp(strcat(['r=',num2str(r),'/',int2str(num_plots+1),', t=',num2str(t),', total dist trav=0, volume=',num2str(liq_volume_temp),', min(R)=',num2str(min(a-Rsol(1,:)))]));
    r=r+1;
    
    % March through time
    while (total_wave_dist<2*pi*la*num_times_after_wave_tracking && t<=tend) 
        num_time_steps = num_time_steps+1;
    
        %Call RHS - 2nd order predictor/corrector method
        Rtemp1 = LW_PRE2012_RHS_S123_Jul2025(R,N,Bo,a,la,da_fact); 
        Rpred=R+dt*Rtemp1;
        
        % Check to see if volume isn't being sufficiently conserved!  
        % If it's not, cut dt by factor of coarse_factor (2)
        % If it is, skip to 'else' statement below
        liq_volume_temp = sum(pi*(a^2-Rpred.^2)*dz);
        mean_R_temp_sq_next = (-liq_volume_temp/(pi*L)+a^2);
        if (abs(mean_R_temp_sq_next-mean_R_temp_sq)>10^-8 && num_dt_halves<num_dt_halves_max) 
            dt=dt/coarse_factor; 
            num_dt_halves=num_dt_halves+1;
            disp('-----------------------------------------------');
            disp(strcat(['dt halving #',int2str(num_dt_halves),', New dt=',num2str(dt),', t=',num2str(t)]));
            disp('-----------------------------------------------');
        else
            Rtemp2 = LW_PRE2012_RHS_S123_Jul2025(Rpred,N,Bo,a,la,da_fact); 
            R = R+dt/2*(Rtemp1+Rtemp2);
    
            %Check volume conservation
            liq_volume_temp = sum(pi*(a^2-R.^2)*dz);
            mean_R_temp_sq = (-liq_volume_temp/(pi*L)+a^2);
    
            % Check how far waves have traveled
            if total_distance_traveled<=2*pi*la*num_times_before_wave_tracking
                [total_distance_traveled temp1_max] = lw_distance_traveled(la, N, R, temp1_max, total_distance_traveled);
            else
            end
    
            % Save snapshots if it's been long enough since the last snapshot
            if t>=time_snaps(r)
                tvec=[tvec t];
                liq_volume = sum(pi*(a^2-R.^2)*dz);
                mean_R(1,r) = (-liq_volume/(pi*L)+a^2)^0.5;
                Rsol(r,:) = a-R(:);
                temporary_minimum=min(a-Rsol(r,:));
                dthalvinglastr=dt;
                temp_crest_loc=[];
                Rtemppad=[R(end-4:end) R R(1:5)];
                for zz=1+5:length(z)+5
                    if Rtemppad(zz)<Rtemppad(zz-1) && Rtemppad(zz)<Rtemppad(zz+1)
                        temp_crest_loc=[temp_crest_loc zz-5];
                    end
                end
                wave_crest_temp_unsort=R(temp_crest_loc);
                [wave_crest_temp idx]=sort(wave_crest_temp_unsort);
                wave_crest_temp_loc=temp_crest_loc(idx);
    
                for qq=1:min([length(wave_crests(1,:)),length(temp_crest_loc)]);
                    wave_crests(r,qq)=R(wave_crest_temp_loc(qq));
                    wave_crest_locs(r,qq)=wave_crest_temp_loc(qq);
                end
    
                padval=round(11/2/dz); %round(pi*sqrt(2)/dz);
                Rtemppad=[R(end-padval+1:end) R R(1:padval)];
                for qq=1:min([length(wave_crests(1,:)),length(temp_crest_loc)]);
                    mean_thicknesses(r,qq)=sqrt(mean((Rtemppad(padval+wave_crest_temp_loc(qq)-padval:padval+wave_crest_temp_loc(qq)+padval)).^2));
                end
    
         
                r = r+1;
            end    
            t = t+dt; 
        end
        if num_dt_halves>=7
            num_refinements=num_refinements+1;
            if num_refinements==1 numtimesrefine='first';
            elseif num_refinements==2 numtimesrefine='second';
            elseif num_refinements==3 numtimesrefine='third';
            elseif num_refinements==4 numtimesrefine='fourth';
            elseif num_refinements==5 numtimesrefine='fifth';
            elseif num_refinements==6 numtimesrefine='sixth';
            elseif num_refinements==7 numtimesrefine='seventh';
            elseif num_refinements==8 numtimesrefine='eighth';
            elseif num_refinements==9 numtimesrefine='ninth';
            end
            disp('***********************************************************');
            disp(strcat(['Refining mesh for ',numtimesrefine,' time.  t=',num2str(t),'; returning to t=',num2str(tvec(end))]));
            disp('***********************************************************');
            %%%%%%%%%% REMOVE THIS BREAK TO INCLUDE MESH REFINEMENT %%%%%%%%%%%
            Rsolmult(:,:,ww)=Rsol;
            wave_crests_mult(:,:,ww)=wave_crests;
            wave_crest_locs_mult(:,:,ww)=wave_crest_locs;
            mean_thicknesses_mult(:,:,ww)=mean_thicknesses;
            
            break
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            t=tvec(end);
            R=a-Rsol(r-1,:);
            dtsave=dtsave/32;
            time_snaps=[time_snaps(1:r-1) time_snaps(r-1)+dtsave:dtsave:time_snaps(r-1)+100*dtsave];
            dt=dthalvinglastr/16;
            num_dt_halves=1+0.1*num_refinements;
            dz=dz/2;
            numptsper2pi=numptsper2pi*2;
            N=N*2;
            z=z(1):dz:z(end)+dz;
            Rfft=fft(R);
            Rfft=[Rfft(1:N/8) zeros(1,3*N/4+1) Rfft(end-N/8+2:end)];
            R=real(ifft(Rfft)*2);
            figure;
            plot(R);
            Rework_Rsol;
            Rsol2=zeros(num_plots+1,N);
            for rrr=1:r-1
                Rfft=fft(Rsol(rrr,:));
                Rfft=[Rfft(1:N/8) zeros(1,3*N/4+1) Rfft(end-N/8+2:end)];
                Rsol2(rrr,:)=real(ifft(Rfft)*2);
            end
            Rsol=Rsol2;
            %Check volume conservation
            liq_volume_temp = sum(pi*(a^2-R.^2)*dz);
            mean_R_temp_sq = (-liq_volume_temp/(pi*L)+a^2);
        end
    end
    
    Rsolmult(:,:,ww)=Rsol;
    wave_crests_mult(:,:,ww)=wave_crests;
    wave_crest_locs_mult(:,:,ww)=wave_crest_locs;
    mean_thicknesses_mult(:,:,ww)=mean_thicknesses;
    
    %progress update to screen
    send(dq, '1');
    
    % record Rsol into Rsol_all
    Rsol_all{ww} = Rsol;

    % record if plugged
    if r<num_plots+1
        plugged(ww) = 1;
        timeToPlug(ww) = t;
    end

    
    disp('Done');

end

% reshape Rsol_all so each row is each individual run, and each column is a different a value
Rsol_all = reshape(Rsol_all, numtests, []);  


%count the numplugs by aggregating plugged vector
numplugs = sum(reshape(plugged, numtests, [])); 
disp(strcat([num2str(numplugs),' runs have ended with a plug']))


fprintf('run time: %f\n', toc);

%%
% Save solution

%get current time
save_time = datestr(now,'yyyymmdd_HHMMSS');

% set data folder name, and save file name
% !!!Remember to change the names accordance to mission number!!!
dataFolderName = 'HPC_test_Data';
fname = sprintf('HPCTest_la12_RsolCollect_%s.mat', save_time);

% Find the data folder that the data is to be saved in
scriptFolder = pwd;                                     %get the full path that the script is in
parentFolder = fileparts(scriptFolder);                 %go up one folder to find the location of the parentfolder
saveFolder = fullfile(parentFolder, dataFolderName);    %set the save folder to be the data folder
if ~exist(saveFolder, 'dir')                            %make a data folder if there is not already one
    mkdir(saveFolder);
end
fnamedotmat= fullfile(saveFolder, [fname,'.mat']);  

%choose what variables to save
save(fnamedotmat, 'plugged','numplugs','timeToPlug','a_values','atests','numtests','la','Bo', ...
    'Rsol_all', 'num_plots',"atests_inc", "atests_max", ...
    'atests_min', 'Ltilde', 'randSeed', '-v7.3');   


%% Tips for processing and graphing

% 1. The Rsol_all stores Rsol from all the runs, each column is trejectories from a
% different 'a' value, and each row is a different trejectory under that same
% 'a' value, but different initial condition

% 2. the a_values is a column vector contains all the a_values that was tested. It is indexed
% the exact same way as the Rsol_all. Ex: Column 3 of Rsol_all all have the
% same 'a' value as a_values(3).

% 3. atests is for the simulation only, an aggregated version of a_value depending on numtest. 
% It is just so the parfor loop can know which 'a' value to use. 
% It's probably for analysis, but we're saving it just in case. 