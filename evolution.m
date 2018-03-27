function [phi1,phi2] = evolution(phi1,phi2,d,epsilon,timestep,mu,v,lamda1,lamda2)
%****************演化函数*********************************


%***************水平集相关函数phi1*******************************
    phi1=NeumannBoundCond(phi1);
    H1 = Heaviside(phi1,epsilon);
    Dirac1 = Dirac(phi1,epsilon);
    Curvature1=curvature(phi1);
    
%*******************水平集相关函数phi2***********************************************

    phi2=NeumannBoundCond(phi2);
    H2 = Heaviside(phi2,epsilon);
    Dirac2 = Dirac(phi2,epsilon);
    Curvature2=curvature(phi2);

%**********************演化方程各项计算***********************************
%      DataForce1=(lamda2*d(:,:,2)-lamda1*d(:,:,1)).*Dirac1;
     lengthterm1=v*Dirac1.*Curvature1;
     dis_regu_term1=mu*Regulation_term(phi1);
     
     lengthterm2=v*Dirac2.*Curvature2;
     dis_regu_term2=mu*Regulation_term(phi2);
     
     DataF1 = -((d(:,:,1)-d(:,:,2)-d(:,:,3)+d(:,:,4)).*H2+d(:,:,2)-d(:,:,4));
     DataF2 = -((d(:,:,1)-d(:,:,2)-d(:,:,3)+d(:,:,4)).*H1+d(:,:,3)-d(:,:,4));
%**********************************************************************

%      phi1 = phi1 +timestep*(DataF.*Dirac1); 
%      phi1 = phi1 +timestep*(DataF.*Dirac1); 

     phi1 = phi1 +timestep*(DataF1.*Dirac1+lengthterm1+dis_regu_term1); 
     phi2 = phi2 +timestep*(DataF2.*Dirac2+lengthterm2+dis_regu_term2); 
     
     
%      phi1 = phi1 +timestep*(DataF1.*Dirac1);%+lengthterm1+dis_regu_term1); 
%      phi2 = phi2 +timestep*(DataF2.*Dirac2);%+lengthterm2+dis_regu_term2); 
     
%     DataF2 = -(d(:,:,1)-d(:,:,2)-d(:,:,3)+d(:,:,4)).*H1-(d(:,:,3)-d(:,:,4));
%     phi2 = phi2 +timestep*(DataF2.*Dirac2);

%************************曲率计算函数**************************
function f=curvature(Phi)
[nx,ny]=gradient(Phi);
s=sqrt(nx.^2+ny.^2);
s=(s==0).*eps+s;
Nx=nx./s;
Ny=ny./s;
[Nxx,Nxy]=gradient(Nx);
[Nyx,Nyy]=gradient(Ny);
f=Nxx+Nyy;
%**************************************************************

%***************距离规则项的计算******************************
function f=Regulation_term(Phi)
[nx,ny]=gradient(Phi);
s=sqrt(nx.^2+ny.^2);
a=(s>=0)&(s<=1);
b=(s>1);
ps=(s-1./((s==0)+(s~=0).*s)).*b+a.*(2*s.^3-3*s.^2+s);
dps=((ps==0)+(ps~=0).*ps)./((s==0)+(s~=0).*s);
[nxx,nxy]=gradient(dps.*nx-nx);
[nyx,nyy]=gradient(dps.*ny-ny);
f=4*del2(Phi)+nxx+nyy;