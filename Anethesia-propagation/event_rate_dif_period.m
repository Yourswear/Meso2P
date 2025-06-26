% For deep and superficial layers during anesthesia induction/recovery, analyze 180s of spike rate data
% Induction: 40-220s (peak synchrony period), Recovery: same window (pre-awakening)
clear; clc;
fs = 7.67;                      % Imaging frame rate (Hz)
timebin = [40 220];             % Analysis time window (seconds)
framebin = floor(timebin*fs);   % Convert to frame indices

%% Load calcium imaging data
dataname = '17_induce_deep';
load(['./',dataname,'.mat']);

%% Parameters
cellnum = size(dff_sig,1);       % Number of neurons
timepoint = size(dff_sig,2);     % Total imaging frames
bin = 9200;                      % Analysis window (2500 or 5000 for other conditions)
sigtoprocess = double(spk_sig(:,1:bin));  % Process spike signals
all_sig = mean(dff_sig(:,1:bin),1);      % Population fluorescence signal
windowsize = 10;                
%% Generate spike time matrix for rate calculation
% Focused analysis: First 2000 frames (subsequent processing requires only 1687)
spike_time_martix = zeros(cellnum,200);
% Create binary spike matrix (1 = spike event)
for i = 1:cellnum
    for j = 1:2000
        spike_time_martix(i,j) = sum(spike_time_estimates{i} == j);
    end
end

% Calculate mean firing rate during analysis window
spk_rate = sum(spike_time_martix(:,framebin(1):framebin(2)), 'all') / ...
        (framebin(2) - framebin(1)) / cellnum * 7.67;  % Convert to spikes/s
disp(spk_rate);  % Display population firing rate