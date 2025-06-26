% For both deep and superficial layers during induce/recover periods, analyze 180s of spike rate data
% Induce: 40-220s (peak synchrony), Recover: same window (pre-awakening)
clear; clc;
fs = 7.67;                      % Imaging frame rate (Hz)
timebin = [40 220];             % Analysis time window (seconds)
framebin = floor(timebin*fs);   % Convert to frame indices
binsize = framebin(2)-framebin(1); % Window size (frames)
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

name = ['pcurve_',period,'_',layer,'_raw.mat'];
datapath = [folderpath,name];
load(datapath);                 % Load pairwise correlation matrix
cellnum = size(pvalue,1);       % Number of neurons

%% Compute mean pairwise correlations 
% Focused analysis: 1687 frames (subset of first 2000 frames)
p_statistic = zeros(cellnum,1);
p_statistic = mean(pvalue(:,framebin(1)+6:23:framebin(2)),2);