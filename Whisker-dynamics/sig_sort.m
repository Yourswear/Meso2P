function sig_sorted = sig_sort(Signal, pointnum)
% Sorts neuronal signals by temporal correlation to stimulus-evoked response template
% Input:
%   Signal:     Neural activity matrix [neurons × timepoints]
%   pointnum:   Estimated number of stimulus events (overridden to 100)
% Output:
%   sig_sorted: Activity matrix sorted by response correlation

pointnum = 100;  % Number of stimulus events for template construction
times = 50;      % Iterations for correlation-based sorting
sz = size(Signal); % [neuron_count, temporal_bins]

%% Preprocessing: Remove silent neurons
for i = 1:sz(1)
    if max(Signal(i,:)) == 0
        Signal(i,:) = [];
        i = i - 1;
    end
end

%% Construct stimulus-aligned template
sumsig = sum(Signal, 1);  % Population activity vector
[B, timespk] = maxk(sumsig, pointnum);  % Identify peak response events
template = zeros(1, sz(2));  % Initialize temporal template
template(timespk) = 255;     % Binarize at stimulus events

%% Initial correlation-based sorting
cormartix = zeros(1, sz(1));  % Correlation coefficients
for i = 1:sz(1)
    cormartix(i) = min(corrcoef(Signal(i,:), template), [], 'all');  % Stimulus-response correlation
end

[C, index] = sort(cormartix);  % Sort by correlation strength
sig_sorted = zeros(sz(1), sz(2));
for i = 1:sz(1)
    sig_sorted(i,:) = Signal(index(i), :);  % Initial sorted matrix
end

%% Iterative template refinement
for i = 1:times
    bina_sig_temp = sig_sorted;  % Current sorted activity
    template_refined = mean(bina_sig_temp(1:10, :), 1);  % Refined template (top 10 neurons)
    
    for j = 1:sz(1)
        cormartix(j) = min(corrcoef(sig_sorted(j,:), template_refined), [], 'all');
    end
    
    [C, index] = sort(cormartix);  % Re-sort by refined correlation
    for j = 1:sz(1)
        sig_sorted(j,:) = bina_sig_temp(index(j), :);  % Update sorting
    end
end
