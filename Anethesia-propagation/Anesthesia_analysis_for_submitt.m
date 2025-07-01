clear; clc; close all;
%% Load data
dataname = '17_induce_deep';
load(['./',dataname,'.mat']);

%% Parameters
cellnum = size(dff_sig,1);       % Number of neurons
timepoint = size(dff_sig,2);     % Total imaging frames
bin = 9200;                      % Analysis window (2500 or 5000 for other conditions)
sigtoprocess = double(spk_sig(:,1:bin));  % Process spike signals
all_sig = mean(dff_sig(:,1:bin),1);      % Population fluorescence signal
windowsize = 10;                 % Spike rate calculate window

%% Population activity analysis
all_spk = zeros(1,bin);          % Initialize population firing rate vector
spike_time_martix = zeros(cellnum,bin);  % Spike time matrix

% Generate binary spike matrix
for i = 1:cellnum
    for j = 1:bin
        spike_time_martix(i,j) = sum(spike_time_estimates{i} == j);
    end
end

% Calculate firing rate (spikes/s)
for i = 1+windowsize/2:bin-windowsize/2
    all_spk(i) = sum(spike_time_martix(:,i-windowsize/2:i+windowsize/2),'all')/...
        (windowsize)/cellnum*7.67; % 7.67 Hz = frame rate
end

% Plot population firing rate
figure(1);
plot(all_spk(1:bin),'LineWidth',2,'color','black'); axis off;
title('Population firing rate');
set(gcf,'Position', [400, 400, 1200, 100]);

%% Neural synchrony analysis
windowsize = 116;  % Correlation window size
step = 23;         % Sliding window step
pcurve = zeros(1,bin);       % Population synchrony curve
pvalue = zeros(cellnum,bin); % Neuron-wise synchrony values

% Calculate pairwise spike correlations
for i = 1+windowsize/2:step:timepoint-windowsize/2
    for j = 1:cellnum
        % Exclude current neuron
        sig_temp = sigtoprocess(:,i-windowsize/2:i+windowsize/2);
        sig_temp(j, :) = [];

        % Compute correlation with population average
        coe_martix = corrcoef(sigtoprocess(j,i-windowsize/2:i+windowsize/2), mean(sig_temp,1));
        pvalue(j,i) = abs(coe_martix(1,2)); % Absolute correlation coefficient
    end
    pcurve(i) = mean(pvalue(:,i),1); % Population synchrony
end

%% Spatiotemporal wave analysis
% Initialize dorsal-ventral activity matrix
tran_hist_deep = zeros(16,bin);
centroid = zeros(cellnum,2);    
for n = 1:cellnum
    centroid(n,:) = properties{n}.centroid;
end
x = centroid(:,1);
y = centroid(:,2);

% Bin neurons by dorsoventral position
for i = 1:16
    tran_hist_deep(i,:) = mean(spk_sig(y>=(i-1)*128 & y<i*128,:),1);
end

% Identify population activity peaks
n = 49; % Example wave index
[pks,locs,w] = findpeaks(mean(spk_sig,1),'Threshold',0.001,'MinPeakHeight',0.04,'MinPeakDistance',20);

% Visualize wave profile
figure(2);
bar3(tran_hist_deep(:,locs(n)-15:locs(n)+15)./max(tran_hist_deep(:,locs(n)-15:locs(n)+15),[],2));
title(['Barplot of spike ',num2str(n)]);
%% Wave propagation visualization
pks_num = length(pks); % Number of detected waves
neuron_deltat_value = zeros(cellnum,1); % Peak latency matrix
image = zeros(2048,2048,pks_num);       % Wave video frames

% Create wave propagation frames
for n = 1:pks_num
    [~,neuron_deltat_value] = max(spk_sig(:,locs(n)-10:locs(n)+10),[],2);
    for i = 1:cellnum
        for j = 1:length(properties{1,i}.coords(:,1))
            % Assign latency values to pixel coordinates
            image(properties{1,i}.coords(j,1)+1,properties{1,i}.coords(j,2)+1,n) = 21 - neuron_deltat_value(i);
        end
    end
end
image(image>20) = 0; % Thresholding
a = movmean(image,5,3); % Temporal smoothing
tiff_save(int16(a),'C:\Users\ElPsyCongroo\Desktop\image.tif'); % Save video

%% Wave propagation dynamics quantification(Global response delay)
time_std_curve = zeros(numel(locs),1); % Propagation variability

% Calculate wave-to-wave propagation consistency
for n = 1:numel(locs)
    % Normalize wave profiles
    tran_hist_norm = tran_hist_deep(1:end-1,locs(n)-10:locs(n)+10)./...
        max(tran_hist_deep(1:end-1,locs(n)-10:locs(n)+10),[],2);
    
    % Find peak latency per dorsoventral bin
    [M, idx] = max(tran_hist_norm,[],2);
    
    % Quantify propagation variability
    time_std_curve(n) = std(idx);
end

% Plot propagation variability over time
figure(3);
plot(locs,mapminmax(movmean(time_std_curve,20)',0,1));
title('Global response delay');
norm_curve = mapminmax(movmean(time_std_curve,20)',0,1)';
t_time = locs/7.67; % Convert frames to seconds