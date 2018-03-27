%*************�����������С��ԭ��minimum relative entropy���Ƶ���������ģ����
%**************����ƫ��У����ָ��ƫ���򿴳��Ǽ��Ե�****************************************************
%******************ʵ������ָ�У��********************************************
%*****************by �����ѧ ɽ��Т    data��2012-01-03*************************
clc;
clear all;
close all;

%*************************����ͼ��******************************
I = imread('095.bmp'); 
I = double(img(:,:,1));
% I = imread('095.bmp');
% I = imread('3Tbrain2.png');

[row,col] = size(I);
% log_I=log(I);
log_I=I;

%**************************************************************

%*******************ת�Ƶ���ֵ��*********************************
% kenel=3;%��ģ���СΪ����1,3,5,7,9............
% edge=(kenel-1);
% pad_num=edge/2;
% 
% img_pad=padarray(log_I,[pad_num pad_num]);%����ģ���С�������
% [row,col]=size(img_pad);
% for i=1:row-edge             %���ѭ���ƶ�ģ��
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
% log_I=img_mean;%�õ�ͼ���ֵ��
% % log_I=uint8(log_I);
% log_I=double(log_I);
% figure(2);
% imagesc(log_I,[0,255]);colormap(gray);
%**************************************************************


%********************����ˮƽ����ʼ��*******************************
% c0=1;
% cw=roipoly;
% phi1=c0*2*(0.5-cw);
% contour(phi1,[0,0],'r');
%******************************************************************

%*********************������˹�˵����******************************
sigma1 =2.5;
K_R=fspecial('gaussian',3,sigma1);%����ˮƽ�������Ĺ��򻯣����ø�˹���򻯶��������ã�

%ģ���ѡȡԭ���ǣ�ǿ���������أ��������أ��Լ��ϳɵ��Ƿ����ʸ߷���ͼ����
%��sigma����ȡС��ͼ��Խ�ǹ⻬��ͬ�ʣ���sigmaȡ�Ľϴ�
sigma2 =8;
K = ones(4*sigma2+1);%����ģ�Ͳ������Ƶĳ�����ģ��
% K=fspecial('gaussian',round(4*sigma2+1),sigma2);%��˹��ģ��
%***********************************************************************


%************˫��ˮƽ����ʼ��**************************************
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

%*****************************�ݻ������ĳ�ʼ��***************8
Iternum =15; %Iterations
epsilon = 1;
timestep = 1;
dim = 4;%����ָ�
lamda1=1;
lamda2=1;
mu=0.1;
v=0.001*255^2;
% v=1;
%**************************************************************************

%*******************%��̬�ֲ��ķ���s��ƫ����b��c�ĳ�ʼ��*****************************
s=ones(row,col,dim).*(1/sqrt(2*pi));
b=ones(row,col);
for i = 1:dim
c(1:row,1:col,i) = i*(max(log_I(:))-min(log_I(:)))./dim + min(log_I(:));
end
%*************************************************************


tic
%*******************************���Ĳ���***************************************
for i = 1:Iternum
    u = compute_u(phi1,epsilon,phi2);% ����heaviside����
    c = compute_c(log_I,K,u,b,s);% ����c 
    b = compute_b(log_I,K,u,c,s);% ����ƫ����b 
    s = compute_s(log_I,b,K,c,u);% ���·���s
    d = computer_d(log_I,K,s,b,c);%���������
    [phi1,phi2]= evolution(phi1,phi2,d,epsilon,timestep,mu,v,lamda1,lamda2);%�����ݻ�����,����Phi
%     phi1 = conv2(phi1,K_R,'same');%ˮƽ���������򻯣����ø�˹���򷽷���
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

%*****************����ʵ�źŽ������********************************
J = zeros(row,col);
for i = 1:dim
    J = J + u(:,:,i).*c(:,:,i);%��ͬ��֯��ǿ����ʵֵ������ʵ�ź�
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

%****************************ƫ������ʾ*****************************
figure(4);
% real_b=exp(double(b));%ָ����ԭ
real_b=b;
imagesc(real_b);colormap(gray);axis off;
title('ƫ����');
%****************************************************************


%****************************��ʵ�ź�(�ָ����)��ʾ*************************
figure(5);
%subplot(1,4,3);    
% J=exp(J);
J=uint8(J);
imagesc(J,[0,255]);colormap(gray);axis off;
title('����ԭʼ�ź�');
%*************************************************************

%*********************ȥ��ƫ������źţ�ƫ��У�������****************************
figure(6);
% corrected_I=exp(log_I-b);
corrected_I=log_I-b;
% b=b+(b==0).*eps;
% corrected_I=uint8(log_I./b);
imagesc(corrected_I);axis off;colormap(gray);
title('ƫ��У����ͼ��');
%********************************************************************

%**********************У��ǰ��ֱ��ͼ�Ա�*********************************
% figure(7);
% subplot(1,2,1);imhist(uint8(I));axis([0 255 0 1000]);title('Hist of Original Image');
% subplot(1,2,2);imhist(uint8(corrected_I));axis([0 255 0 1000]);title('Hist of MRE');colormap(gray);
%**************************************************************************

%*************У��ǰ��ֱ��ͼ�Ա�**********************
% figure(8);
% [counts1,x1]=imhist(uint8(I));[counts2,x2]=imhist(uint8(corrected_I));
% subplot(1,2,1);plot(x1,counts1);axis([0 255 0 1400]);title('Hist of Original Image');colormap(gray);
% subplot(1,2,2);plot(x2,counts2);axis([0 255 0 1400]);title('Hist of MREVLS');colormap(gray);
%**************************************************

%***************************GM��GM����ʾ***************************************
img1=u(:,:,1).*corrected_I;
img2=u(:,:,2).*corrected_I;%u2��GM
img3=u(:,:,3).*corrected_I;
img4=u(:,:,4).*corrected_I;%u4��WG
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

%*******************************У�����GM��WM�۵��Ա�******************************
% figure(12);
% [count_img1,x_img1]=imhist(uint8(img1));[count_img3,x_img3]=imhist(uint8(img33));
% subplot(2,2,1);plot(x_img1,count_img1,'r');axis([0 255 0 700]);title('Hist of corrected WM');colormap(gray);
% subplot(2,2,3);plot(x_img3,count_img3,'b');axis([0 255 0 700]);title('Hist of corrected GM');colormap(gray);
% subplot(2,2,[2,4]);plot(x_img1,count_img1,'r',x_img3,count_img3,'b');axis([0 255 0 700]);title('corrected Overlap between GM and WM');colormap(gray);
%****************************************************************

%*****************����У�����CJVϵ��************************************
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
% CJV_corrected=(WM_stdVar+GM_stdVar)/abs(WM_mean-GM_mean)        %���CJVϵ����


%*******************У��ǰ��GM��WM�۵��Ա�******************************************
% img1_orginal=u(:,:,1).*log_I;
% img2_orginal=u(:,:,2).*log_I;%u2��GM
% img3_orginal=u(:,:,3).*log_I;
% img4_orginal=u(:,:,4).*log_I;%u4��WG
% img33_orginal=img2_orginal+img3_orginal;
% figure(13);
% [count_img1_orginal,x_img1_orginal]=imhist(uint8(img1_orginal));[count_img3_orginal,x_img3_orginal]=imhist(uint8(img33_orginal));
% subplot(2,2,1);plot(x_img1_orginal,count_img1_orginal,'r');axis([0 255 0 500]);title('Hist of Original WM');colormap(gray);
% subplot(2,2,3);plot(x_img3_orginal,count_img3_orginal,'b');axis([0 255 0 500]);title('Hist of Original GM');colormap(gray);
% subplot(2,2,[2,4]);plot(x_img1_orginal,count_img1_orginal,'r',x_img3_orginal,count_img3_orginal,'b');axis([0 255 0 500]);title('Original Overlap between GM and WM');colormap(gray);
%*************************************************************************

%*********************����У��ǰ��CJVϵ��***************************
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
% CJV_orignal=(WM_stdVar+GM_stdVar)/abs(WM_mean-GM_mean)        %���CJVϵ����
%**************************************************************************










