clc;
clear;
%%
folderpath='D:\Private-HuJiaHao\for_calcium_data_process\file_cai-1';
cell = 'cell1'
index='_001';
period = 11500;
subindex=['_00',num2str(period)];
%加载GT电生理记录
load([folderpath,'\processed_data\data_20120521_',cell,index,'.mat']);
%加载新型流程处理结果
load([folderpath,'\Spike_inferrence_test\',cell,index,subindex,'\',cell,index,subindex,'.mat']);
% %加载OASIS处理结果
load([folderpath,'\Spike_inferrence_test\',cell,index,subindex,'\',cell,index,subindex,'_oasis_sig.mat']);
%电生理记录结果导入
t_ephys=obj.timeSeriesArrayHash.value{4}.time;
detected_spikes=obj.timeSeriesArrayHash.value{5}.valueMatrix;
spike_time=t_ephys(detected_spikes);
%bin spikes
spike_binned=histcounts(spike_time,[0,obj.timeSeriesArrayHash.value{1}.time']);
spike_binned = spike_binned((period-1)*2400+1:period*2400);
% if subindex=='_001'
%     spike_binned=spike_binned(1:2400);
%     spike_binned=spike_binned(2401:4800);
% end
%smooth bin spikes
standard_deviation_smoothing = 2.4; % in number of frames; can be chosen depending on what is considered optimal or desired
spike_binned_smoothed = conv(spike_binned,fspecial('gaussian',[40 1],standard_deviation_smoothing),'same');
% cas_sig_binned = histcounts(spike_time_estimates,[1:2401]);
% cas_sig_smoothed = conv(cas_sig_binned,fspecial('gaussian',[40 1],standard_deviation_smoothing),'same');

% disp(corrcoef(OASISSignal,spike_binned_smoothed(1,1:2400)));
disp(corrcoef(spk_sig(33:2367),spike_binned_smoothed(1,33:2367)));
disp(corrcoef(spk_oasis(33:2367),spike_binned_smoothed(1,33:2367)));
figure(777);
plot(spike_binned_smoothed);
hold on
plot(spk_sig);
% hold on
% plot(OASISSignal./30);