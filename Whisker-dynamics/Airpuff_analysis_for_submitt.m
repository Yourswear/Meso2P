clear; clc; close all;
%% Load data
dataname = '17-AH-all_result';
load(['./',dataname,'.mat']);
load(['./','AP0,6mm_recode','.mat']);

%% Parameters
cellnum = size(dff_sig,1);
timepoint = size(dff_sig,2);
AP = 0;  % Anterior-Posterior coordinate index (0 for bregma)
bin = [61 90];  % Time window for analysis (frames)
binsize = bin(2)-bin(1)+1;
bia = 462;  % Inter-trial interval (frames)

%% Part 1: Brain area distribution analysis
% Create neuron mask from coordinates
mask_neuron = zeros(2048,2048);
for i = 1:cellnum
    for j = 1:length(properties{i}.coords)
        mask_neuron(properties{i}.coords(j,1)+1,properties{i}.coords(j,2)+1) = 255;
    end
end
mask = fliplr(L');
mask = circshift(mask,40);
mask = circshift(mask,10,2);
preview = imadd(mask*8,uint8(mask_neuron)); figure(1); imshow(preview);

% Define brain areas based on AP coordinate
if AP == -1
    areanum = 16;
    areaname = {'ACA' 'MOs','MOp','SSp_ul','SSp_II','SSp_un','SSp_bfd',...
    'SSp_tr','VISa','VISam','VISpm','VISrl','VISp','RSPagl','RSPd','RSPv'};
elseif AP == 0  % Bregma position
    areanum = 13;
    areaname = {'PL','ACA','MOs','MOp','SSp_ul','SSp_II','SSp_un',...
        'SSp_tr','VISa','VISam','RSPagl','RSPd','RSPv'};
    M_index = [3,4];      % Motor cortex indices
    S_index = [5,6,7,8];  % Somatosensory cortex indices
    V_index = [9,10];     % Visual cortex indices
    RSP_index = [11,12,13]; % Retrosplenial cortex indices
end

% Identify active cells based on firing rate
active_cell = [];
for i = 1:length(spike_time_estimates)
    baseline = sum(spike_time_estimates{i}>0&spike_time_estimates{i}<2000)/1999;
    active_level = 0;
    for j = 1:10
            active_level = active_level + sum((spike_time_estimates{i}>2298+(j-1)*462+1)& ...
            (spike_time_estimates{i}<2298+(j-1)*462+115));
    end
    active_level = active_level/1150;
    if active_level>2*baseline
        active_cell = [active_cell,i];
    end
end

%% Create variables for each brain area and count neurons
for i = 1:areanum
    eval([areaname{i},'=','[];']);
end
for i = 1:cellnum
    whicharea = mask(floor(properties{i}.centroid(1)),floor(properties{i}.centroid(2)));
    if whicharea~=0
        eval([areaname{whicharea},'=','[',areaname{whicharea},',',num2str(i),'];']);
    end
end

%% Store area distribution in cell array
area_distribution = cell(2,areanum);
for i = 1:areanum
    area_distribution{1,i} = areaname{i};
    area_distribution{2,i} = eval(areaname{i});
end

%% Plot neuron count per area
neurons = zeros(1,areanum);
for i = 1:areanum
    neurons(i) = numel(area_distribution{2,i});
    area_distribution{3,i} =  neurons(i);
end
figure(2); bar(categorical(area_distribution(1,:)), neurons);
area_distribution = area_distribution';

%% Part 2: Consolidate all trial signals into single matrix for heatmap
time1 = zeros(10,1); time2 = zeros(10,1);
cellnum = length(spike_time_estimates);
spk_alltrial = zeros(10*cellnum,268);  % Preallocate trial matrix

% Extract airpuff response period
for trial = 1:10
    time1(trial) = (trial-1)*bia + 2298 - 77 + 1;
    time2(trial) = (trial-1)*bia + 2298 + 115 + 76;
    spk_alltrial((trial-1)*cellnum+1:trial*cellnum,:) = spk_sig(:,time1(trial):time2(trial));
end
figure(3), imagesc(spk_alltrial); clim([-0.05 0.5]);

%% Part 3: Baseline vs trial spike rates (bar plot, units: spk/s)
time1 = time1 + bin(1);
time2 = time1 + bin(2) - bin(1);

% Create spike time matrix
spike_time_martix = zeros(cellnum,timepoint);
for i = 1:cellnum
    for j = 1:timepoint
        spike_time_martix(i,j) = sum(spike_time_estimates{i} == j);
    end
end

% Calculate spike rates during airpuff trials
spk_rate_10trials = zeros(10,1);
for trial = 1:10
    spk_rate_10trials(trial) = sum(spike_time_martix(:,time1(trial):time2(trial)),'all')/...
        (binsize)/cellnum*7.67; % 7.67 = imaging frame rate
end

% Calculate baseline spike rate
spk_rate_base = sum(spike_time_martix(:,1:2000),'all')/2000/cellnum*7.67;
spk_rate_all = [spk_rate_base; spk_rate_10trials];
figure(4); bar(spk_rate_all);

%% Part 4: Temporal activity heatmap across trials
tran_hist_all_trial_temp = zeros(16,30,10);
centroid = zeros(cellnum,2);    
for n = 1:cellnum
    centroid(n,:) = properties{n}.centroid;
end
x = centroid(:,1);
y = centroid(:,2);
for trial = 1:10
    for i = 1:16
        tran_hist_all_trial_temp(i,:,trial) = mean(spk_sig(y>=(i-1)*128&y<i*128,time1(trial):time2(trial)),1);
    end
end

% Align to stimulus onset
tran_hist_all_trial = zeros(16,15,10);
for i = 1:10
    [M, idx] = max(tran_hist_all_trial_temp(:,15:23,i),[],2);
    center = mode(idx+14);
    tran_hist_all_trial(:,:,i) = tran_hist_all_trial_temp(:,center-7:center+7,i);
end
tran_hist_all_trial = tran_hist_all_trial./max(tran_hist_all_trial,[],2);

%% Plot temporal activity for trial 1
figure(5); subplot(2,1,1);
set(gcf,'Position', [700, 0, 700, 1400]);
x = linspace(3,-3,16)';  % Dorsal-ventral axis
y = linspace(-0.75,0.75,30)'; % Time axis
barfig = bar3(x,tran_hist_all_trial(:,:,1)); view(220,50);
set(barfig,'LineWidth',2);
set(gca,'Gridcolor','black'); set(gca,'GridLineWidth',0.5');
set(gca, 'GridAlpha', 1);
set(gca, 'Fontname','Arial','FontSize',16);
ylim([-3.5 3.5]);
hx = xlabel('', 'FontSize', 16, 'FontName', 'Arial','Rotation',43);
hy = ylabel('', 'FontSize', 16, 'FontName', 'Arial','Rotation', -45);
set(hx, 'Position', [5, -5, -2]);
set(hy, 'Position', [-5, 5, -2]);
set(gca, 'LineWidth', 2);
set(gca, 'XTickLabel', {'-1', '0', '1'});
set(gca, 'XTick', linspace(1, 15, 3)); 
set(gca, 'YTick', linspace(-3, 3, 3)); 
colormap sky;
daspect([2,1,0.15]);

%% Plot temporal activity for trial 10
subplot(2,1,2);
x = linspace(3,-3,16)';
y = linspace(-0.75,0.75,30)';
barfig = bar3(x,tran_hist_all_trial(:,:,10)); view(220,50);
set(barfig,'LineWidth',2);
set(gca,'Gridcolor','black'); set(gca,'GridLineWidth',0.5');
set(gca, 'GridAlpha', 1);
set(gca, 'Fontname','Arial','FontSize',16);
ylim([-3.5 3.5]);
hx = xlabel('', 'FontSize', 16, 'FontName', 'Arial','Rotation',43);
hy = ylabel('', 'FontSize', 16, 'FontName', 'Arial','Rotation', -45);
set(hx, 'Position', [5, -5, -2]);
set(hy, 'Position', [-5, 5, -2]);
set(gca, 'LineWidth', 2);
set(gca, 'XTickLabel', {'-1', '0', '1'});
set(gca, 'XTick', linspace(1, 15, 3)); 
set(gca, 'YTick', linspace(-3, 3, 3)); 
colormap sky;
daspect([2,1,0.15]);

%% Part 5: Trial-to-trial pattern r using CCA
valid = [1,2,3,4,5,6,8,9,10];  % Valid trial indices
bin = [73,84]; % Time window for CCA (±1.5s around stimulus)

% Define airpuff response period
for trial = 1:10
    time1(trial) = (trial-1)*bia + 2298 - 77 + 1;
    time2(trial) = (trial-1)*bia + 2298 + 115 + 76;
end
binsize = bin(2) - bin(1) + 1;
time1 = time1 + bin(1);
time2 = time1 + bin(2) - bin(1);

% Compute CCA r matrix
CCA_trial_martix = zeros(10,10); % (trial x trial similarity)
for trial1 = 1:10
    for trial2 = 1:10
        X = spk_sig(:,time1(trial1):time2(trial1));
        Y = spk_sig(:,time1(trial2):time2(trial2));
        [wxMat,wyMat,rVec] = SparseCCA(X',Y',1,1,2,1);
        CCA_trial_martix(trial1,trial2) = rVec;
    end
end
temp = [];
for i = 1:10
    for j = 1:10
        if i ~= j
            temp = [temp, CCA_trial_martix(i,j)];
        end
    end
end
maxvalue = max(temp,[],'all');
minvalue = min(temp,[],'all');
for i = 1:10
    for j = 1:10
        if i ~= j
            CCA_trial_martix(i,j) = (CCA_trial_martix(i,j) - minvalue)/(maxvalue - minvalue);
        end
    end
end

%% Part 5b: Brain area-specific pattern 
% Extract region-specific signals
M = []; S = []; V = []; RSP = [];
for i = 1:length(M_index)
    M = [M; spk_sig(area_distribution{M_index(i),2},:)];
end
for i = 1:length(S_index)
    S = [S; spk_sig(area_distribution{S_index(i),2},:)];
end
for i = 1:length(V_index)
    V = [V; spk_sig(area_distribution{V_index(i),2},:)];
end
for i = 1:length(RSP_index)
    RSP = [RSP; spk_sig(area_distribution{RSP_index(i),2},:)];
end
cellnum_M = size(M,1); cellnum_S = size(S,1);
cellnum_V = size(V,1); cellnum_RSP = size(RSP,1);

% Create region-specific trial matrices
spk_alltrial_M = zeros(10*cellnum_M,binsize);
spk_alltrial_S = zeros(10*cellnum_S,binsize);
spk_alltrial_V = zeros(10*cellnum_V,binsize);
spk_alltrial_RSP = zeros(10*cellnum_RSP,binsize);

for trial = 1:10
    spk_alltrial_M((trial-1)*cellnum_M+1:trial*cellnum_M,:) = M(:,time1(trial):time2(trial));
end
for trial = 1:10
    spk_alltrial_S((trial-1)*cellnum_S+1:trial*cellnum_S,:) = S(:,time1(trial):time2(trial));
end
for trial = 1:10
    spk_alltrial_V((trial-1)*cellnum_V+1:trial*cellnum_V,:) = V(:,time1(trial):time2(trial));
end
for trial = 1:10
    spk_alltrial_RSP((trial-1)*cellnum_RSP+1:trial*cellnum_RSP,:) = RSP(:,time1(trial):time2(trial));
end

figure('Name','M'), imagesc(spk_alltrial_M); clim([0 0.5]);
figure('Name','S'), imagesc(spk_alltrial_S); clim([0 0.5]);
figure('Name','V'), imagesc(spk_alltrial_V); clim([0 0.5]);
figure('Name','RSP'), imagesc(spk_alltrial_RSP); clim([0 0.5]);

% Compute CCA r for each region
CCA_trial_martix_M = zeros(10,10,5);
CCA_trial_martix_S = zeros(10,10,5);
CCA_trial_martix_V = zeros(10,10,5);
CCA_trial_martix_RSP = zeros(10,10,5);

for trial1 = 1:10
    for trial2 = 1:10
        X = spk_alltrial_M((trial1-1)*cellnum_M+1:trial1*cellnum_M,:);
        Y = spk_alltrial_M((trial2-1)*cellnum_M+1:trial2*cellnum_M,:);
        [wxMat,wyMat,rVec] = SparseCCA(X',Y',1,1,2,1);
        CCA_trial_martix_M(trial1,trial2,:) = rVec;
    end
end

for trial1 = 1:10
    for trial2 = 1:10
        X = spk_alltrial_S((trial1-1)*cellnum_S+1:trial1*cellnum_S,:);
        Y = spk_alltrial_S((trial2-1)*cellnum_S+1:trial2*cellnum_S,:);
        [wxMat,wyMat,rVec] = SparseCCA(X',Y',1,1,2,1);
        CCA_trial_martix_S(trial1,trial2,:) = rVec;
    end
end

for trial1 = 1:10
    for trial2 = 1:10
        X = spk_alltrial_V((trial1-1)*cellnum_V+1:trial1*cellnum_V,:);
        Y = spk_alltrial_V((trial2-1)*cellnum_V+1:trial2*cellnum_V,:);
        [wxMat,wyMat,rVec] = SparseCCA(X',Y',1,1,2,1);
        CCA_trial_martix_V(trial1,trial2,:) = rVec;
    end
end

for trial1 = 1:10
    for trial2 = 1:10
        X = spk_alltrial_RSP((trial1-1)*cellnum_RSP+1:trial1*cellnum_RSP,:);
        Y = spk_alltrial_RSP((trial2-1)*cellnum_RSP+1:trial2*cellnum_RSP,:);
        [wxMat,wyMat,rVec] = SparseCCA(X',Y',1,1,2,1);
        CCA_trial_martix_RSP(trial1,trial2,:) = rVec;
    end
end
figure(6); set(gcf,'Position', [0, 400, 2000, 200]);
subplot(1,5,1); heatmap(CCA_trial_martix); clim([0.4 0.8]); title('All');
subplot(1,5,2); heatmap(CCA_trial_martix_M(:,:,1)); clim([0.4 0.8]); title('Motor');
subplot(1,5,3); heatmap(CCA_trial_martix_S(:,:,1)); clim([0.4 0.8]); title('Somatosensory');
subplot(1,5,4); heatmap(CCA_trial_martix_V(:,:,1)); clim([0.4 0.8]); title('Visual');
subplot(1,5,5); heatmap(CCA_trial_martix_RSP(:,:,1)); colormap jet; clim([0.4 0.8]); title('Retrosplenial');

%% Part 6: Cross-region CCA projection analysis
% Initialize change matrices
S_changemartix = zeros(cellnum_S,10);
M_changemartix = zeros(cellnum_M,10);
V_changemartix = zeros(cellnum_V,10);
RSP_changemartix = zeros(cellnum_RSP,10);

for which_trial = 1:10
    % Extract trial data for each region
    S_trial = S(:,time1(which_trial):time2(which_trial));
    M_trial = M(:,time1(which_trial):time2(which_trial));
    V_trial = V(:,time1(which_trial):time2(which_trial));
    RSP_trial = RSP(:,time1(which_trial):time2(which_trial));
    
    % Compute cross-region CCA
    [VSx,VSy,rVec] = SparseCCA(V_trial',S_trial',1,1,2,1);
    [VMx,VMy,rVec] = SparseCCA(V_trial',M_trial',1,1,2,1);
    [VRSPx,VRSPy,rVec] = SparseCCA(V_trial',RSP_trial',1,1,2,1);
    [SMx,SMy,rVec] = SparseCCA(M_trial',S_trial',1,1,2,1);
    [SRSPx,SRSPy,rVec] = SparseCCA(RSP_trial',S_trial',1,1,2,1);
    [MRSPx,MRSPy,rVec] = SparseCCA(M_trial',RSP_trial',1,1,2,1);
    
    for mode = 1  % Use first CCA mode
        % Compute region-specific projections
        S_CCA_sig = (S_trial'*VSy(:,mode)*VSy(:,mode)' + ...
                     S_trial'*SMy(:,mode)*SMy(:,mode)' + ...
                     S_trial'*SRSPy(:,mode)*SRSPy(:,mode)')'./3;
        M_CCA_sig = (M_trial'*VMy(:,mode)*VMy(:,mode)' + ...
                     M_trial'*SMx(:,mode)*SMx(:,mode)' + ...
                     M_trial'*MRSPx(:,mode)*MRSPx(:,mode)')'./3;
        V_CCA_sig = (V_trial'*VSx(:,mode)*VSx(:,mode)' + ...
                     V_trial'*VMx(:,mode)*VMx(:,mode)' + ...
                     V_trial'*VRSPx(:,mode)*VRSPx(:,mode)')'./3;
        RSP_CCA_sig = (RSP_trial'*VRSPy(:,mode)*VRSPy(:,mode)' + ...
                       RSP_trial'*SRSPx(:,mode)*SRSPx(:,mode)' + ...
                       RSP_trial'*MRSPy(:,mode)*MRSPy(:,mode)')'./3;
        
        % Store CCA weights
        S_changemartix(:,which_trial) = (abs(VSy) + abs(SMy) + abs(SRSPy))./3;
        M_changemartix(:,which_trial) = (abs(SMx) + abs(VMy) + abs(MRSPx))./3;
        V_changemartix(:,which_trial) = (abs(VSx) + abs(VMx) + abs(VRSPx))./3;
        RSP_changemartix(:,which_trial) = (abs(SRSPx) + abs(MRSPy) + abs(VRSPy))./3;
        
        % Plot region-specific CCA projections
        figure(7); subplot(10,1,which_trial); imagesc(S_CCA_sig); box off; axis off; clim([0 0.5])
        figure(8); subplot(10,1,which_trial); imagesc(M_CCA_sig); box off; axis off; clim([0 0.5])
        figure(9); subplot(10,1,which_trial); imagesc(V_CCA_sig); box off; axis off; clim([0 0.5])
        figure(10); subplot(10,1,which_trial); imagesc(RSP_CCA_sig); box off; axis off; clim([0 0.5])
    end
end
S_changemartix = sig_sort(S_changemartix,6);
M_changemartix = sig_sort(M_changemartix,6);
V_changemartix = sig_sort(V_changemartix,6);
RSP_changemartix = sig_sort(RSP_changemartix,6);

%% Part 7: PCA analysis of cross-trial changes
figure(11);
subplot(1,4,1); imagesc(M_changemartix); title('Motor');
subplot(1,4,2); imagesc(S_changemartix); title('Somatosensory');
subplot(1,4,3); imagesc(V_changemartix); title('Visual');
subplot(1,4,4); imagesc(RSP_changemartix); title('Retrosplenial');

% Perform PCA on cross-trial patterns
[coeff, scoreS, latent, tsquared, explained, mu] = pca(S_changemartix');
[coeff, scoreM, latent, tsquared, explained, mu] = pca(M_changemartix');
[coeff, scoreV, latent, tsquared, explained, mu] = pca(V_changemartix');
[coeff, scoreRSP, latent, tsquared, explained, mu] = pca(RSP_changemartix');

% Plot first principal component across trials
figure(12);
subplot(4,1,1); plot(scoreS(valid,1)); title('Somatosensory PC1');
subplot(4,1,2); plot(scoreM(valid,1)); title('Motor PC1');
subplot(4,1,3); plot(scoreV(valid,1)); title('Visual PC1');
subplot(4,1,4); plot(scoreRSP(valid,1)); title('Retrosplenial PC1');