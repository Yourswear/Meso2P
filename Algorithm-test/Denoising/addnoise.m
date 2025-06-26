%该函数用于给raw data添加高斯-泊松混合噪声
clc;
clear;
sample_fre = [1,2,3,4,5,10,15,20,25,30];%[5,10,15,20,25,30];
folderpath='D:\Lab-Share\Member-HuJiaHao\20240704-GUI测试数据\2-不同帧率\不同帧率原始数据\512\';
savepath='D:\Lab-Share\Member-HuJiaHao\20240704-GUI测试数据\2-不同帧率\不同帧率加了噪声之后的数据\512\';
for i = sample_fre
    Image = double(tiffreadVolume([folderpath,num2str(i),'_Hz_512','.tif']));
    sz=size(Image);
    sigma=100;%高斯噪声的标准差参数
    for SNR=100 %通过改变Q值来改变每个像素值对应的光子数，例如像素值为10，实际收集光子数可能是1000也可能是5000
        tic
        %Q=10^(SNR/10)/260;
        Q=0.001
        for j=1:sz(3)
            pos = poissrnd(Q.*Image(:,:,j));
            pos=pos./Q;
            gau_pos=pos+normrnd(0,sigma,sz(1),sz(2));
            noised=uint16(gau_pos); 
            imwrite(noised,[savepath,'0dB_',num2str(i),'_Hz_512_noised','.tif'],'WriteMode','append');
        end
        toc
    end
end


%Image=uint16(Image);
%write(Image)
