%绘图程序，把GT和不同算法(SUNS,Suite2p,CaImAn)输出的结果放在同一张图上,cellpose的绘图在代码Comp_cellpose中
clc;
clear;
%% gt加载 输出gt
gtpath = 'D:\Lab-Share\Member-HuJiaHao\大视场项目\20250505_GUIcomp\';
gtname = 'GT_mask_properties.mat';
load([gtpath,gtname]);
gt_pro = properties;
gt = zeros(290,290);
for i = 1:length(gt_pro)
    for j = 1:length(gt_pro{i}.coords)
        gt(gt_pro{i}.coords(j,1)+1,gt_pro{i}.coords(j,2)+1)=1;
    end
end
figure(1);
imshow(gt);
%% NeuroPixelAI处理结果加载 输出nai
naipath = 'D:\Lab-Share\Member-HuJiaHao\大视场项目\20250505_GUIcomp\不同信噪比分割对比\masks_suite2p\';
nainame = '18.4dB';
load([naipath,nainame,'.mat']);
sui_pro = stat;
sui = zeros(290,290);
for i = 1:length(sui_pro)
    for j = 1:length(sui_pro{i}.xpix)
        sui(sui_pro{i}.ypix(j)+1,sui_pro{i}.xpix(j)+1)=1;
    end
end
figure(2);
imshow(sui);

comp_mask = zeros(290,290,3);
comp_mask(:,:,1)=50*gt;
comp_mask(:,:,2)=50*sui;
imshow(comp_mask);
