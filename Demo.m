%*************本程序采用最小熵原理（minimum relative entropy）推导出的能量模型来
%**************进行偏差校正与分割，将偏差域看成是加性的****************************************************
%******************实现四相分割校正********************************************
%*****************by山金孝 *************************
clc;
clear all;
close all;

%*************************输入图像******************************
I = imread('095.bmp'); 
I = double(img(:,:,1));
% I = imread('095.bmp');
% I = imread('3Tbrain2.png');

[row,col] = size(I);
% log_I=log(I);
log_I=I;

%**************************************************************

%*******************转移到均值域*********************************
% kenel=3;%核模板大小为奇数1,3,5,7,9............
% edge=(kenel-1);
% pad_num=edge/2;
% 
% img_pad=padarray(log_I,[pad_num pad_num]);%根据模板大小扩充矩阵
% [row,col]=size(img_pad);
% for i=1:row-edge             %外层循环移动模板
%     for j=1:col-edge
%         
%         padarray_sub=img_pad([i:i+edge],[j:j+edge]);
%         padarray_sub_mean=sum(sum(padarray_sub))/kenel^2;
%         
%         img_mean(i,j)=padarray_sub_mean;
%     end
% end
% kenel=5;
% k_mean=ones(kenel);
% a=conv2(I,k_mean,'same');
% img_mean=a./kenel^2;
% log_I=img_mean;%得到图像均值域
% % log_I=uint8(log_I);
% log_I=double(log_I);
% figure(2);
% imagesc(log_I,[0,255]);colormap(gray);
%**************************************************************


%********************单相水平集初始化*******************************
% c0=1;
% cw=roipoly;
% phi1=c0*2*(0.5-cw);
% contour(phi1,[0,0],'r');
%******************************************************************

%*********************两个高斯核的设计******************************
sigma1 =2.5;
K_R=fspecial('gaussian',3,sigma1);%用以水平集函数的规则化，采用高斯规则化抖动大，少用！

%模板的选取原则是：强度异质严重，噪声严重（自己合成的那幅异质高方差图），
%则sigma尽量取小；图像越是光滑和同质，则sigma取的较大
sigma2 =8;
K = ones(4*sigma2+1);%用以模型参数估计的常量核模板
% K=fspecial('gaussian',round(4*sigma2+1),sigma2);%高斯核模板
%***********************************************************************


%************双相水平集初始化**************************************
% phi1=ones(row,col);
% phi2=ones(row,col);
% 
% o_y1=round(row/2);o_x1=round(col/2)-round(col/8);
% o_y2=round(row/2);o_x2=round(col/2)+round(col/8);
% r1=round(row/8);r2=round(row/8);
% 
% for i=1:row
%     for j=1:col
%         if ((i-o_y1)^2+(j-o_x1)^2)<r1^2
%             phi1(i,j)=-1;
%         end
%         if ((i-o_y2)^2+(j-o_x2)^2)<r2^2
%             phi2(i,j)=-1;
%         end
%     end
% end

load biasButNoNoise_cross.mat
% load phi_cross_3t.mat
phi1=phi_1;
phi2=phi_2;

figure(1);
imagesc(I,[0,255]);colormap(gray);hold on;axis off;
contour(phi1,[0 0],'r','LineWidth',2);
contour(phi2,[0 0],'b','LineWidth',2);
%*********************************************************************

%*****************************演化参数的初始化***************8
Iternum =15; %Iterations
epsilon = 1;
timestep = 1;
dim = 4;%四相分割
lamda1=1;
lamda2=1;
mu=0.1;
v=0.001*255^2;
% v=1;
%**************************************************************************

%*******************%正态分布的方差s，偏差域b和c的初始化*****************************
s=ones(row,col,dim).*(1/sqrt(2*pi));
b=ones(row,col);
for i = 1:dim
c(1:row,1:col,i) = i*(max(log_I(:))-min(log_I(:)))./dim + min(log_I(:));
end
%*************************************************************


tic
%*******************************核心部分***************************************
for i = 1:Iternum
    u = compute_u(phi1,epsilon,phi2);% 更新heaviside函数
    c = compute_c(log_I,K,u,b,s);% 更新c 
    b = compute_b(log_I,K,u,c,s);% 更新偏差域b 
    s = compute_s(log_I,b,K,c,u);% 更新方差s
    d = computer_d(log_I,K,s,b,c);%数据项计算
    [phi1,phi2]= evolution(phi1,phi2,d,epsilon,timestep,mu,v,lamda1,lamda2);%调用演化方程,更新Phi
%     phi1 = conv2(phi1,K_R,'same');%水平集函数规则化（采用高斯规则方法）
%     phi2 = conv2(phi1,K_R,'same');
    if(mod(i,1)==0)
        figure(2);
        pause(0.0001);
    
        imagesc(I,[0 255]);colormap(gray);axis off;
        hold on;
        contour(phi1,[0 0],'r');%,'LineWidth',2);
%         contour(phi1,[0 0],'r','LineWidth',2);
    if dim == 4
       contour(phi2,[0 0],'b');%,'LineWidth',2);
%        contour(phi2,[0 0],'b','LineWidth',2);
    end   
       iterNum=[num2str(i), ' iterations'];  
       title(iterNum);hold off;   
    end
end
%*********************************************************************
toc

%*****************对真实信号进行组合********************************
J = zeros(row,col);
for i = 1:dim
    J = J + u(:,:,i).*c(:,:,i);%不同组织的强度真实值，即真实信号
end
%*******************************************************************

%-------------------------
% figure(3);
% %subplot(1,4,1);
% imagesc(I,[0 255]);axis off;colormap(gray);
% hold on;
% contour(phi1,[0 0],'r');
% % contour(phi2,[0 0],'b');
% hold off;
%subplot(1,4,2);

%****************************偏差域显示*****************************
figure(4);
% real_b=exp(double(b));%指数复原
real_b=b;
imagesc(real_b);colormap(gray);axis off;
title('偏差域');
%****************************************************************


%****************************真实信号(分割后结果)显示*************************
figure(5);
%subplot(1,4,3);    
% J=exp(J);
J=uint8(J);
imagesc(J,[0,255]);colormap(gray);axis off;
title('区域原始信号');
%*************************************************************

%*********************去除偏差域后信号（偏差校正结果）****************************
figure(6);
% corrected_I=exp(log_I-b);
corrected_I=log_I-b;
% b=b+(b==0).*eps;
% corrected_I=uint8(log_I./b);
imagesc(corrected_I);axis off;colormap(gray);
title('偏差校正后图像');
%********************************************************************

%**********************校正前后直方图对比*********************************
% figure(7);
% subplot(1,2,1);imhist(uint8(I));axis([0 255 0 1000]);title('Hist of Original Image');
% subplot(1,2,2);imhist(uint8(corrected_I));axis([0 255 0 1000]);title('Hist of MRE');colormap(gray);
%**************************************************************************

%*************校正前后直方图对比**********************
% figure(8);
% [counts1,x1]=imhist(uint8(I));[counts2,x2]=imhist(uint8(corrected_I));
% subplot(1,2,1);plot(x1,counts1);axis([0 255 0 1400]);title('Hist of Original Image');colormap(gray);
% subplot(1,2,2);plot(x2,counts2);axis([0 255 0 1400]);title('Hist of MREVLS');colormap(gray);
%**************************************************

%***************************GM与GM的显示***************************************
img1=u(:,:,1).*corrected_I;
img2=u(:,:,2).*corrected_I;%u2是GM
img3=u(:,:,3).*corrected_I;
img4=u(:,:,4).*corrected_I;%u4是WG
% figure(9);
% imagesc(img1,[0,255]);colormap(gray);axis off;
% figure(10);
% imagesc(img2,[0,255]);colormap(gray);axis off;
img33=img3+img2;
figure(10);
imagesc(img33,[0,255]);colormap(gray);axis off;
figure(11);
imagesc(img4,[0,255]);colormap(gray);axis off;
%***********************************************************************

%*******************************校正后的GM与WM折叠对比******************************
% figure(12);
% [count_img1,x_img1]=imhist(uint8(img1));[count_img3,x_img3]=imhist(uint8(img33));
% subplot(2,2,1);plot(x_img1,count_img1,'r');axis([0 255 0 700]);title('Hist of corrected WM');colormap(gray);
% subplot(2,2,3);plot(x_img3,count_img3,'b');axis([0 255 0 700]);title('Hist of corrected GM');colormap(gray);
% subplot(2,2,[2,4]);plot(x_img1,count_img1,'r',x_img3,count_img3,'b');axis([0 255 0 700]);title('corrected Overlap between GM and WM');colormap(gray);
%****************************************************************

%*****************计算校正后的CJV系数************************************
% GM_var=0;
% GM_mean=0;
% WM_var=0;
% WM_mean=0;
% img_WM=(img1>50).*img1;img_GM=(img33>50).*img33;
% N=0;M=0;s1=0;s2=0;
% for i=1:row
%     for j=1:col
%         if(img_WM(i,j)~=0)
%             N=N+1;
%         end
%         if(img_GM(i,j)~=0)
%             M=M+1;
%         end
%     end
% end
% WM_mean=sum(sum(img_WM))/N;
% GM_mean=sum(sum(img_GM))/M;
% 
% for i=1:row
%     for j=1:col
%         if(img_WM(i,j)~=0)
%             s1=s1+(img_WM(i,j)-WM_mean)^2;
%         end
%         if (img_GM(i,j)~=0)
%             s2=s2+(img_GM(i,j)-GM_mean)^2;
%         end
%     end
% end
% WM_stdVar=sqrt(s1/N);
% GM_stdVar=sqrt(s2/M);
% 
% CJV_corrected=(WM_stdVar+GM_stdVar)/abs(WM_mean-GM_mean)        %输出CJV系数；


%*******************校正前的GM与WM折叠对比******************************************
% img1_orginal=u(:,:,1).*log_I;
% img2_orginal=u(:,:,2).*log_I;%u2是GM
% img3_orginal=u(:,:,3).*log_I;
% img4_orginal=u(:,:,4).*log_I;%u4是WG
% img33_orginal=img2_orginal+img3_orginal;
% figure(13);
% [count_img1_orginal,x_img1_orginal]=imhist(uint8(img1_orginal));[count_img3_orginal,x_img3_orginal]=imhist(uint8(img33_orginal));
% subplot(2,2,1);plot(x_img1_orginal,count_img1_orginal,'r');axis([0 255 0 500]);title('Hist of Original WM');colormap(gray);
% subplot(2,2,3);plot(x_img3_orginal,count_img3_orginal,'b');axis([0 255 0 500]);title('Hist of Original GM');colormap(gray);
% subplot(2,2,[2,4]);plot(x_img1_orginal,count_img1_orginal,'r',x_img3_orginal,count_img3_orginal,'b');axis([0 255 0 500]);title('Original Overlap between GM and WM');colormap(gray);
%*************************************************************************

%*********************计算校正前的CJV系数***************************
% GM_var=0;
% GM_mean=0;
% WM_var=0;
% WM_mean=0;
% img_WM=(img1_orginal>50).*img1_orginal;img_GM=(img33_orginal>50).*img33_orginal;
% N=0;M=0;s1=0;s2=0;
% for i=1:row
%     for j=1:col
%         if(img_WM(i,j)~=0)
%             N=N+1;
%         end
%         if(img_GM(i,j)~=0)
%             M=M+1;
%         end
%     end
% end
% WM_mean=sum(sum(img_WM))/N;
% GM_mean=sum(sum(img_GM))/M;
% 
% for i=1:row
%     for j=1:col
%         if(img_WM(i,j)~=0)
%             s1=s1+(img_WM(i,j)-WM_mean)^2;
%         end
%         if (img_GM(i,j)~=0)
%             s2=s2+(img_GM(i,j)-GM_mean)^2;
%         end
%     end
% end
% WM_stdVar=sqrt(s1/N);
% GM_stdVar=sqrt(s2/M);
% 
% CJV_orignal=(WM_stdVar+GM_stdVar)/abs(WM_mean-GM_mean)        %输出CJV系数；
%**************************************************************************










