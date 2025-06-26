clc;
clear;
noisedpath='D:\Lab-Share\Member-HuJiaHao\20240704-GUI测试数据\2-不同帧率\不同帧率加噪声_deepinter\\';
cleanpath='D:\Lab-Share\Member-HuJiaHao\20240704-GUI测试数据\2-不同帧率\不同帧率原始数据\512\\';
for i=[1,2,3,4,5,10,15,20,25,30]
    cleanname=[num2str(i),'_Hz_512'];%cleanname=[num2str(i),'_Hz'];
    clean290 = double(tiffreadVolume([cleanpath,cleanname,'.tif']));
    clean290=clean290(:,:,31:130);
    %clean290 = rescale(double(tiffreadVolume([cleanpath,cleanname,'.tif'])),0,5000);
    %clean288=clean290(2:end-1,2:end-1,:);
    sigPower=sum(clean290.^2,'all');%信号总能量计算
    noised = double(tiffreadVolume([noisedpath,'0dB_',num2str(i),'Hz_512_deepinter_output','.tif']));
    noised=noised(:,:,1:100);
    noisePower=sum(abs(clean290-noised).^2,'all');%噪声总能量计算
    Real_SNR=10*log10(sigPower/noisePower);
    disp([num2str(Real_SNR)]);
end

