% Given the multiple fluctuations in intracellular calcium concentration 
% during anesthesia, a sliding time window was employed to calculate ΔF/F0 
% for long-term signal processing to ensure analytical accuracy, 
% based on specific characteristics of the data.
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
    
%% Calculate dF/F
win = 300;  % Baseline window size
dff = zeros(cellnum,timepoint);

% Compute dF/F using percentile baseline
for i = 1+win/2:timepoint-win/2
    bas = raw_sig(:,i-win/2:i+win/2);
    bas_sort = sort(bas,2);  % Sort fluorescence values
    baseline = mean(bas_sort(:,1:5),2);  % Use lowest 5% as baseline
    dff(:,i) = ((double(raw_sig(:,i)) - baseline)) ./ baseline;
end
dff_sig = dff;
