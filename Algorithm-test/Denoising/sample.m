clear;
clc;

sample_fre = 3:5;%[5,10,15,20,25,30];
folderpath = 'D:\Lab-Share\Member-HuJiaHao\test\';
savepath = 'D:\Lab-Share\Member-HuJiaHao\20240704-GUI测试数据\2-不同帧率\';

for i = sample_fre
    sample_interval = 100./i;
    sample_Image = uint16(zeros(290,290,floor(30000/(100/i))));
    for j = 1:9
        Image = 1000.*double(tiffreadVolume([folderpath,'100Hz-20240706_0000',num2str(j),'.tif']));
        for k = (j-1)*(2500/sample_interval)+(1:2500/sample_interval)
            sample_Image(:,:,k) = Image(:,:,floor((k-(j-1)*2500/sample_interval)*sample_interval));
            imwrite(sample_Image(:,:,k),[savepath,num2str(i),'_Hz.tif'],'WriteMode','append');
        end
    end
    for j = 10:12
        Image = 1000.*double(tiffreadVolume([folderpath,'100Hz-20240706_000',num2str(j),'.tif']));
        for k = (j-1)*(2500/sample_interval)+(1:2500/sample_interval)
            sample_Image(:,:,k) = Image(:,:,floor((k-(j-1)*2500/sample_interval)*sample_interval));
            imwrite(sample_Image(:,:,k),[savepath,num2str(i),'_Hz.tif'],'WriteMode','append');
        end
    end
end
