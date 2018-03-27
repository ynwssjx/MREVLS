%************************7T头颅水平集初始化**************************
clc;
clear all;
close all;
row=256;col=256;
data=ones(row,col);
phi_1=data;
phi_2=data;
for y=1:row
    for x=1:col
        if (y>=50)&(y<=190)
            s1=(y-50)*(-100/170)+150;
            s2=(y-50)*(-80/-170)+150;
            if(x>=s1)&(x<=s2)
                phi_1(y,x)=-1;
            end
        end
    end
end
% contour(phi_1,[0,0],'r');

for y=1:row
    for x=1:col
        if (y>=50)&(y<=280)
            w1=(y-50)*(80/170)+80;
            w2=(y-50)*(-80/170)+210;
            if(x>=w1)&(x<=w2)
                phi_2(y,x)=-1;
            end
        end
    end
end
contour(phi_1,[0,0],'r');hold on;
contour(phi_2,[0,0],'r');

% save phi_1 phi_2 7T_cross;load
save('biasButNoNoise_cross','phi_1','phi_2')








