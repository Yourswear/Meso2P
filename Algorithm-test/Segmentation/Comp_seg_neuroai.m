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
naipath = 'D:\Lab-Share\Member-HuJiaHao\大视场项目\20250505_GUIcomp\masks_neuroai\';
nainame = '-6dB_Enhanced_result';
load([naipath,nainame,'.mat']);
nai_pro = properties;
nai = zeros(290,290);
for i = 1:length(nai_pro)
    for j = 1:length(nai_pro{i}.coords)
        nai(nai_pro{i}.coords(j,1)+1,nai_pro{i}.coords(j,2)+1)=1;
    end
end
figure(2);
imshow(nai);

comp_mask = zeros(290,290,3);
comp_mask(:,:,1)=50*gt;
comp_mask(:,:,2)=50*nai;
imshow(comp_mask);
