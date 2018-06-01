clc
clear
format long
set(0,'FormatSpacing','compact')
Node([1 2],:)=[1 -488.0021  729.2153;2 -507.7515  729.2153] ;
xp=[-645.2600 693.0400]';
ElNo=1;
NE(ElNo,:)=[0     1];
NU=0.2;
SG=1e9/2/(1+NU);
Hi=zeros(2,4);
LE(ElNo)=19.7494;
Den=2200.0000 ;
AG=[0 -9.8]';
 SV11=0
 SV12=0;
 SV22=0;
% 高斯求解
GaussPt=5;
 syms x;         %计算高斯点坐标及系数
maple('with','orthopoly');
f=maple('P',GaussPt,x);		%n为阶数
c=sym2poly(f);
x = sort(roots(c));
p = legendre(GaussPt-1,x);
c = 2*(1-x.^2)./(GaussPt*p(1,:)').^2;
GaussXW = [x c]';
Gi=zeros(2,4);
for iG=1:GaussPt
    x=GaussXW(1,iG);
    xl=([1-x x+1]/2*Node([1 2],2:end)).';
    r=xl-xp;
    rl=sqrt(r.'*r)  ;        
    drd1=r(1)/rl;
    drd2=r(2)/rl;
    drdn=[drd1 drd2]*NE(ElNo,:)';
    NShape=[(1-x)/2 0 (1+x)/2 0;0 (1-x)/2 0 (1+x)/2];
   SV11=SV11+1/8/pi*(LE(ElNo)/2)*GaussXW(2,iG)*(1/(1-NU)*(NU*(2*NE(ElNo,:)*[drd1 drd2]'*Den*AG'*[drd1 drd2]'+...
         (1-2*log(1/rl))*Den*AG'*NE(ElNo,:)')-Den*AG'*[drd1 drd2]'*(2*NE(ElNo,1)*drd1)));
   SV12=SV12+1/8/pi*(LE(ElNo)/2)*GaussXW(2,iG)*(2*NE(ElNo,:)*[drd1 drd2]'*(Den*AG(2)*drd1)+...
        1/(1-NU)*(-Den*AG'*[drd1 drd2]'*(NE(ElNo,:)*[drd2 drd1]')+...
         (1-2*NU)/2*(1-2*log(1/rl))*(Den*AG(2)*NE(ElNo,1))));
   SV22=SV22+1/8/pi*(LE(ElNo)/2)*GaussXW(2,iG)*(2*NE(ElNo,:)*[drd1 drd2]'*(2*Den*AG(2)*drd2)+1/(1-NU)*(NU*((2*NE(ElNo,:)*[drd1 drd2]')*Den*AG'*[drd1 drd2]'+...
        (1-2*log(1/rl))*(Den*AG'*NE(ElNo,:)'))-Den*AG'*[drd1 drd2]'*(2*NE(ElNo,2)*drd2)+(1-2*NU)/2*(1-2*log(1/rl))*(2*Den*AG(2)*NE(ElNo,2))));
end
 [SV11 SV12 SV22]
%解析解
               SV11=0;
                 SV12=0;
                 SV22=0;
                b=NE(ElNo,1);
                Angle=acos(NE(ElNo,1));
                c=cos(pi-Angle);
                if NE(ElNo,2)>=0
                    a=cos(Angle-pi/2);
                else
                    a=cos(Angle+pi/2);
                end
%                     建立一个局部坐标的数列
                angleR=acos(NE(ElNo,:)*[0 1]');
                if NE(ElNo,1)>=0 
                    Trans=[cos(angleR) -sin(angleR);sin(angleR) cos(angleR)];
                else
                    Trans=[cos(angleR) sin(angleR);-sin(angleR) cos(angleR)];
                end 
                LocalNode=zeros(size(Node));
                LocalNode(1,:)=[1 (Trans*(Node(1,2:3)'-xp))'];
                LocalNode(2,:)=[2 (Trans*(Node(2,2:3)'-xp))'];
                
                x2=LocalNode(2,2);
                x1=LocalNode(1,2);
                y0=LocalNode(2,3);
                n1=NE(1);
                n2=NE(2);
                    CP=LE(ElNo)/8/pi/(x2-x1);
                    BGrav1=(x2-x1)-y0*(atan(x2/y0)-atan(x1/y0));
                    BGrav2=(log(x2^2+y0^2)-log(x1^2+y0^2))*y0/2;
                    BGrav3=y0*(atan(x2/y0)-atan(x1/y0));
                    SC1=a*c*BGrav1+(a^2+b*c)*BGrav2+a*b*BGrav3;
                    SC2=c^2*BGrav1+2*a*c*BGrav2+a^2*BGrav3;
                    SC3=a^2*BGrav1+2*a*b*BGrav2+b^2*BGrav3;
                    TC1=(x2-x1)+x2*log(x2^2+y0^2)-x1*log(x1^2+y0^2)-2*BGrav1;
                    S11=1/(1-NU)*(NU*(2*n1*Den*AG(2)*SC1+2*n2*Den*AG(2)*SC2+n2*Den*AG(2)*TC1)-2*n1*Den*AG(2)*SC1);
                    S12=2*Den*AG(2)*(n1*SC3+n2*SC1)+1/(1-NU)*(-Den*AG(2)*(n1*SC2+n2*SC1)+(1-2*NU)/2*Den*AG(2)*n1*TC1);
                    S22=4*Den*AG(2)*(n1*SC1+n2*SC2)+1/(1-NU)*(NU*(2*Den*AG(2)*(n1*SC1+n2*SC2)+n2*Den*AG(2)*TC1)...
                        -2*n2*Den*AG(2)*SC2+(1-2*NU)*Den*AG(2)*n2*TC1);
                    SV11=S11*CP;
                    SV12=S12*CP;
                    SV22=S22*CP;
                    [SV11 SV12 SV22]
%    for iG=1:GaussPt
%        xi=GaussXW(1,iG);
%        x=(1-xi)/2*x1+(1+xi)/2*x2;
%        tyq=tyq+x^3/(x^2+y0^2)^3*GaussXW(2,iG);
%    end
%    tyq=tyq*(x2-x1)/2

