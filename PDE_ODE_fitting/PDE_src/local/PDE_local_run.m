%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reed Ogrosky, Aya Yu
% VCU
% September 2019
% Updated September 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Rsol, time_snaps, tvec, z] = PDE_sim(R, tend, num_plots, la, Bo, a)
    close all;
    
    % Set spatial resolution

    L = 2*pi*la;              % Don't change - Domain is la*(2*pi)
    numptsper2pi = 32;        % Probably don't need to change - sets spatial resolution - number of z gridpoints per 2pi
    N = numptsper2pi*la;      % Don't change - Total number of z gridpoints over entire domain
    dz = L/N;                 % Don't change - Space between successive z gridpoints
    z = 0:dz:L-dz;            % Don't change - Vector containing all z gridpoints
    num_refinements=0;        % Don't change for now; only reason to increase is to see free surface closer to plug formation
    
    % Set initial temporal resolution 
    t0 = 0;                                  % Don't change - Start time
    tnumpersec1=ceil(16000*(a-1)^3);         % Number of time steps per time unit (initially)
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
        Rtemp1 = LW_PRE2012_RHS_S123_Jul2025(R,N,Bo,a,la, da_fact); 
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
                % for zz=1+5:length(z)+5
                %     if Rtemppad(zz)<Rtemppad(zz-1) && Rtemppad(zz)<Rtemppad(zz+1)
                %         temp_crest_loc=[temp_crest_loc zz-5];
                %     end
                % end
                % wave_crest_temp_unsort=R(temp_crest_loc);
                % [wave_crest_temp idx]=sort(wave_crest_temp_unsort);
                % wave_crest_temp_loc=temp_crest_loc(idx);
                % 
                % for qq=1:min([length(wave_crests(1,:)),length(temp_crest_loc)]);
                %     wave_crests(r,qq)=R(wave_crest_temp_loc(qq));
                %     wave_crest_locs(r,qq)=wave_crest_temp_loc(qq);
                % end
    
                % padval=round(11/2/dz); %round(pi*sqrt(2)/dz);
                % Rtemppad=[R(end-padval+1:end) R R(1:padval)];
                % for qq=1:min([length(wave_crests(1,:)),length(temp_crest_loc)]);
                %     mean_thicknesses(r,qq)=sqrt(mean((Rtemppad(padval+wave_crest_temp_loc(qq)-padval:padval+wave_crest_temp_loc(qq)+padval)).^2));
                % end
    
                % Print progress update to screen
                disp(strcat(['r=',num2str(r),'/',int2str(num_plots+1),', t=',num2str(t),', trav=',num2str(total_distance_traveled),', volume=',num2str(liq_volume_temp),', mean(R)=',num2str(mean_R_temp_sq),', min(R)=',num2str(temporary_minimum)]));
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
            % Check volume conservation
            liq_volume_temp = sum(pi*(a^2-R.^2)*dz);
            mean_R_temp_sq = (-liq_volume_temp/(pi*L)+a^2);
        end
    end
    disp('Done');
end

%% Generate a movie from the snapshots saved
function play_movie(Rsol, a, la, num_plots)
    L = 2*pi*la;              % Don't change - Domain is la*(2*pi)
    numptsper2pi = 32;        % Probably don't need to change - sets spatial resolution - number of z gridpoints per 2pi
    N = numptsper2pi*la;      % Don't change - Total number of z gridpoints over entire domain
    dz = L/N;                 % Don't change - Space between successive z gridpoints
    z = 0:dz:L-dz;        



    ymin = 0.0;                           % Smallest value of h to show
    ymax = a;                             % Largest value of h to show
    fnamemov=strcat('Si1'); 
    v=VideoWriter(fnamemov,'MPEG-4');
    v.FrameRate=20; % Set frame rate
    open(v);        % Once you run open(v), you won't be able to change the properties of the video (like framerate, etc.)
    fig1 = figure;
    winsize = get(fig1,'Position'); 
    winsize(1:2) = [0 0]; 
    last_row=0;
    for n=1:1:num_plots+1;
        if (Rsol(n,1)==0 && last_row==0)
            last_row=n-1;
        end
    end
    if last_row>0
        num_frames=last_row;
    else
        num_frames=num_plots+1;
    end
    eval([strcat(sprintf('Movie')) ' = moviein(num_frames,fig1,winsize);']);
    set(fig1,'NextPlot','replacechildren');
    
    % Plot each snapshot
    for r=1:1:num_frames
        plot(z,Rsol(r,:),'k');
        hold on;
        plot(z,2*a-Rsol(r,:),'k');
        fill([z(1) z z(end)],[0 Rsol(r,:) 0],[0.5 0.5 1]);
        fill([z(1) z z(end)],[2*a 2*a-Rsol(r,:) 2*a],[0.5 0.5 1]);
        hold off
        xlabel('z');
        ylabel('h=a-R');
        set(gca,'YTick',[0,a,2*a]);
        set(gca,'YTickLabels',{'-a','0','a'})
        axis([0,L,2*ymin,2*ymax]);
    
        frame=getframe(gcf);
        %writeVideo(v,frame);
    end
    for r=1:50
        hold on;
        s4=fill([z(1) z(1) z(end) z(end)],[2*a a a 2*a],[1 1 1],'Linestyle','none');
        s4.FaceAlpha=0.05; 
        frame=getframe(gcf);
        %writeVideo(v,frame);
    end
    %save movie_interface;
    close(v);
end

%% Parameters
% Model Parameters
a=1.34;                                  % Film thickness ratio (tube radius to core radius) - must be real number greater than 1
Bo=1;                                    % Bond number
la= 12;                                    % Enter a positive integer to set domain size

% Simulation Parameters
tend = 2200; %50*round(1/(a-1)^3);       % End time
num_plots =400;                          % Number of solution snapshots to save 

%% Construct initial condition - R=1+some small perturbations

L = 2*pi*la;              % Don't change - Domain is la*(2*pi)
numptsper2pi = 32;        % Probably don't need to change - sets spatial resolution - number of z gridpoints per 2pi
N = numptsper2pi*la;      % Don't change - Total number of z gridpoints over entire domain
dz = L/N;                 % Don't change - Space between successive z gridpoints
z = 0:dz:L-dz;            % Don't change - Vector containing all z gridpoints

beta=0.01;               % amplitude of small perturbations
num_modes= 5*la;     % Number of Fourier modes to add to initial condition
R=1;                     % Don't change - Mean value of R
for m=1:1:num_modes
    R=R+(a-1)*beta*cos(m*z/la+10000*rand*pi/2);
end

%% Run
[Rsol, time_snaps, tvec] = PDE_sim(R, tend, num_plots, la, Bo, a);

%% Run multiple times
% num_sims = 10;  % number of times of simulation to run
% Rsol_tens = zeros(num_plots+1, (32*la), num_sims);              % initialize 3d tensor to store each simulation in a differnt layer
% for i = 1:num_sims
%     % Construct initial condition - R=1+some small perturbations
%     L = 2*pi*la;              % Don't change - Domain is la*(2*pi)
%     numptsper2pi = 32;        % Probably don't need to change - sets spatial resolution - number of z gridpoints per 2pi
%     N = numptsper2pi*la;      % Don't change - Total number of z gridpoints over entire domain
%     dz = L/N;                 % Don't change - Space between successive z gridpoints
%     z = 0:dz:L-dz;            % Don't change - Vector containing all z gridpoints
% 
%     beta=0.01;               % amplitude of small perturbations
%     num_modes= 1; %5*la;     % Number of Fourier modes to add to initial condition
%     R=1;                     % Don't change - Mean value of R
%     for m=1:1:num_modes
%         R=R+(a-1)*beta*cos(m*z/la+10000*rand*pi/2);
%     end
% 
%     %Run
%     [Rsol, time_snaps, tvec] = PDE_sim(R, tend, num_plots, la, Bo, a);
%     Rsol_tens(:, :, i) = Rsol;
% end
%Rsol_avg = mean(Rsol_tens, 3);


%% Find Rmin, Rmin_ind, Calculate speed of wave
Rmtrx = -Rsol + a; % save all the R values from Rsol
[Rmin Rmin_ind] = min(Rmtrx, [], 2);  % Find the minimum at each time snapshot and save both R and index
Rspeed = zeros(length(Rmin_ind) - 1, 1);
time_inc = time_snaps(2) - time_snaps(1);  % time increment
for i = 2:1:length(Rmin)
    d1 = (Rmin_ind(i) - Rmin_ind(i-1)) * dz;
    d2 = (Rmin_ind(i)*dz + L) - (Rmin_ind(i-1) * dz);
    dtrue = min([abs(d1), abs(d2)]);
    Rspeed(i-1) = dtrue/time_inc;
end

%% Find relative mins


for i = 1:size(Rsol_all, 2)
    Rsol = Rsol_all{i};
    Rmtrx = -1 .* Rsol + 1.35; 
    figure;
    hold on;
    for row = 1:size(Rmtrx, 1)
        for col = 2:size(Rmtrx, 2)-1
            if Rmtrx(row, col) < Rmtrx(row, col - 1) && Rmtrx(row, col) < Rmtrx(row, col + 1)
                scatter(row, Rmtrx(row, col), '.');
            end
        end
    
    end
    hold off;
end

%% Save solution
dataFolderName = 'Flow_fit_Data';
timestamp = datestr(now,'yyyy-mm-dd_HH-MM-SS');
fname = ['PDE_fit_test' timestamp];
% Find the data folder that the data is to be saved in
scriptFolder = pwd;        %get the full path that the script is in
parentFolder = fileparts(scriptFolder);                 %go up one folder to find the location of the parentfolder
saveFolder = fullfile(parentFolder, dataFolderName);    %set the save folder to be the data folder
if ~exist(saveFolder, 'dir')                            %make a data folder if there is not already one
    mkdir(saveFolder);
end

fnamedotmat= fullfile(saveFolder, [fname,'.mat']);                     
save(fnamedotmat,'la','Bo', 'a', 'Rsol', 'Rmin_ind', 'tvec', 'time_snaps');    %choose what variables to save


%% Plot the crest height
figure;
plot(time_snaps(1:end-1), Rmin(1:end-1) );

%% Plot Crest Speed vs Height

figure;
plot(Rmin(1:end-1), Rspeed(1:end));
xlabel('R_{min}');
ylabel('Wave Horizontal Speed');

%% Plot average crest height from multi-r
% figure;
% hold on;
% for ii = 1:num_sims
%     Rsol = Rsol_tens(:, :, ii);
% 
%     Rmtrx = -Rsol + a; % save all the R values from Rsol
%     Rmin = min(Rmtrx, [], 2);  % Find the minimum at each time snapshot, that is the crest
%     plot(time_snaps(1:end-1), Rmin(1:end-1) );
% end