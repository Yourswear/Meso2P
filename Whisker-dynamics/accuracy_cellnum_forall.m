% Input:
%   A: ncells x ntrials1 x nframes
%   B: ncells x ntrials2 x nframes
% 
% Output:
%   p: Decoding accuracy
% 
clear; clc; close all;

%% Data information
folderpath = 'D:\Lab-Share\Member-HuJiaHao\大视场项目\气吹\20250319吹正脸\17\result_20250413\';
name = '17-AH-all_result';

%% Data load---spike_time_estimates
load([folderpath,name,'.mat']);
sz = size(spk_sig);
spk = spk_sig;
neur_num = size(spk,1);

% Initialize accuracy matrix for different cell counts
accu_cellnum = zeros(10,10);
t = 0;

% Loop over different numbers of neurons (from 200 to 4000 in steps of 200)
for ncell = 200:200:4000%floor(0.05*neur_num):floor(0.05*neur_num):neur_num
    disp(['t=',num2str(t)]);
    t = t+1;
    
    % Repeat each cell count 10 times with random neuron selection
    for n = 1:10
        disp(['n=',num2str(n)]);
        % Randomly select ncell neurons
        nspk = spk(randperm(neur_num, ncell),:);
        
        %% Select and preprocess data
        % Input:
        %   A: ncells x ntrials1 x nframes
        %   B: ncells x ntrials2 x nframes
        ntrials1 = 9;
        ntrials2 = 10;
        ncells = size(nspk,1);
        
        % Initialize data arrays
        A = zeros(ncells,ntrials1,200);
        B = zeros(ncells,ntrials2,200);
        
        % Extract trial data for condition A
        for i = 1:ntrials1
            A(:,i,:) = nspk(:,(i-1)*460+2298-85+1:(i-1)*460+2298+115);
        end
        
        % Extract trial data for condition B
        for i = 1:ntrials2
            B(:,i,:) = nspk(:,(i-1)*200+501:(i-1)*200+700);
        end

        %% Decoding analysis
        binsize = 16;       % Temporal bin size for feature extraction
        nsamples = 16;      % Number of bootstrap samples
        nframes = size(A,3)-binsize;  % Total number of time frames
        
        % Initialize accuracy matrix
        accu = nan(nsamples,nframes);
        
        % Main decoding loop
        for sample = 1:nsamples
            for frame = 1:nframes
                % Extract time-binned data
                x = A(:,:,frame:frame+binsize);
                y = B(:,:,frame:frame+binsize);
                
                % Compute decoding accuracy using SVM
                accu(sample,frame) = svm_decoder(x,y,0); 
            end
        end
        
        % Calculate mean accuracy during stimulus period
        accu_cellnum(n,t) = mean(accu(:,85-binsize/2:85+115-binsize),'all');
        
        % Plot temporal accuracy profile
        plot(mean(accu,1)); 
    end
end

% Plot final accuracy vs. number of neurons
plot(mean(accu_cellnum,1));