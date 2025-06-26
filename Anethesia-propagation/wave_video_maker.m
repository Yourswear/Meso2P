clear; clc;
dataname = '17_induce_deep';
load(['./',dataname,'.mat']);

%% Parameters
cellnum = size(dff_sig,1);       % Number of neurons
timepoint = size(dff_sig,2);     % Total imaging frames
bin = 9200;                      % Analysis window (2500 or 5000 for other conditions)
sigtoprocess = double(spk_sig(:,1:bin));  % Process spike signals
all_sig = mean(dff_sig(:,1:bin),1);      % Population fluorescence signal
windowsize = 10;                 % Spike rate calculate window

%% Generate population activity trace
all_spk = zeros(1,bin);          % Population firing rate vector
spike_time_martix = zeros(cellnum,bin);  % Binary spike matrix

% Create spike time matrix (1 = spike event)
for i = 1:cellnum
    for j = 1:bin
        spike_time_martix(i,j) = sum(spike_time_estimates{i} == j);
    end
end

% Calculate smoothed firing rate (spikes/s)
for i = 1+windowsize/2:bin-windowsize/2
    all_spk(i) = sum(spike_time_martix(:,i-windowsize/2:i+windowsize/2),'all')/...
        (windowsize)/cellnum*7.67; % 7.67 Hz = frame rate
end

% Identify population activity peaks
figure('Name','Population firing rate');
plot(all_spk);
findpeaks(all_spk,'Threshold',0.001,'MinPeakHeight',0.06,'MinPeakDistance',20);
[pks,locs,w] = findpeaks(mean(spk_sig,1),'Threshold',0.001,'MinPeakHeight',0.04,'MinPeakDistance',20);

%% Generate wave propagation frames
pks_num = length(pks);           % Number of detected waves
neuron_deltat_value = zeros(cellnum,1); % Peak latency matrix
image = zeros(2048,2048,pks_num);       % Wave video frames

% Create wave propagation frames
for n = 1:pks_num
    % Determine per-neuron peak latency within ±10 frames of population peak
    [~,neuron_deltat_value] = max(spk_sig(:,locs(n)-10:locs(n)+10),[],2);
    
    % Map latency values to neuron coordinates
    for i = 1:cellnum
        for j = 1:length(properties{1,i}.coords(:,1))
            image(properties{1,i}.coords(j,1)+1,properties{1,i}.coords(j,2)+1,n) = ...
                21 - neuron_deltat_value(i); % Invert for visualization
        end
    end
end

% Post-processing and video generation
image(image>20) = 0;             % Remove out-of-range values
a = movmean(image,5,3);          % Temporal smoothing (5-frame window)
tiff_save(int16(a),'C:\Users\ElPsyCongroo\Desktop\image.tif'); % Save as TIFF stack