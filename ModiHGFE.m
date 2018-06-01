function [Hi Gi VFi]=ModiHGFE(GaussXW,Node,LEC,Mat,AG)
% ENCo,xp,GaussXW,GaussPt
    np=Node(1,1);	En1=Node(2,1);  En2=Node(3,1);    xp=Node(1,2:3)';    
	Den=Mat(1);    NU=Mat(3);    SG=Mat(4);
    Hi=zeros(2,4);	Gi=zeros(2,4);	VFi=zeros(2,1);
    
     DLE=Node(3,2:3)-Node(2,2:3);
     LEE=sum(DLE.*DLE,2).^0.5;
     DLE=DLE./[LEE LEE];
     NEE=DLE*[0 -1;1 0];
     NE=NEE;

%  需要进行局部坐标转换 Localxp=[0 0];
%                 定义的常数
% b=NE(ElNo,1);
b=NE(1);
Angle=acos(NE(1));
c=cos(pi-Angle);
if NE(2)>=0
    a=cos(Angle-pi/2);
else
    a=cos(Angle+pi/2);
end
%                     建立一个局部坐标的数列
angleR=acos(NE*[0 1]');
if NE(1)>=0 
    Trans=[cos(angleR) -sin(angleR);sin(angleR) cos(angleR)];
else
    Trans=[cos(angleR) sin(angleR);-sin(angleR) cos(angleR)];
end 
LocalNode=zeros(size(Node));
LocalNode(2,:)=[1 (Trans*(Node(2,2:3)'-xp))'];
LocalNode(3,:)=[2 (Trans*(Node(3,2:3)'-xp))'];

x2=LocalNode(3,2);
x1=LocalNode(2,2);
y0=LocalNode(3,3);
% 参数
n1=NE(1);
n2=NE(2);
% 源点在单元上
%   Gi的求解
if En1==np |En2==np
    AX=diff(Node(2:3,2))/2;
    AY=diff(Node(2:3,3))/2;
    SR=(AX^2+AY^2)^0.5;  
    DE=4*pi*SG*(1-NU);               
    SR2=2*SR;   %一个单元的长度
    DE2=2*DE;
    XXU=3-4*NU;
    X21=2*AX;
    Y21=2*AY;
    %   有重力时需用伽辽金张量表示的带刚体位移的基本解
    temp1=SR*XXU*(1.5-log(SR2))-(7-8*NU)/2*SR;
    temp2=SR*XXU*(0.5-log(SR2))-(7-8*NU)/2*SR;
    AG1=X21^2/(2*SR2)+temp1;
    AG4=Y21^2/(2*SR2)+temp1;
    BG1=X21^2/(2*SR2)+temp2;
    BG4=Y21^2/(2*SR2)+temp2;
    AG2=X21*Y21/(2*SR2);
    if En2==np
        temp=AG1;   AG1=BG1;    BG1=temp;
        temp=AG4;   AG4=BG4;    BG4=temp;                       
    end 
    Gi =[AG1 AG2 BG1 AG2;AG2 AG4 AG2 BG4]/DE2; 
else
    Constp=LEE/8/pi/SG/(1-NU)/(x2-x1)^2;
    temp1=-(7-8*NU)/2*(x2*(x2-x1)-(x2^2-x1^2)/2);
    temp2=-(7-8*NU)/2*(-x1*(x2-x1)+(x2^2-x1^2)/2);
    Treat2=x2*((x2-x1)-y0*(atan(x2/y0)-atan(x1/y0)))...
         -((x2^2-x1^2)-y0^2*log((x2^2+y0^2)/(x1^2+y0^2)))/2;
    Treat6=-x1*((x2-x1)-y0*(atan(x2/y0)-atan(x1/y0)))...
              +((x2^2-x1^2)/2-y0^2/2*log((x2^2+y0^2)/(x1^2+y0^2)));
    y2Const2=x2*y0*(atan(x2/y0)-atan(x1/y0))-y0^2*(log(x2^2+y0^2)-log(x1^2+y0^2))/2;
    y2Const7=-x1*y0*(atan(x2/y0)-atan(x1/y0))+y0^2*(log(x2^2+y0^2)-log(x1^2+y0^2))/2;
    yConst6=y0*x2/2*log((x2^2+y0^2)/(x1^2+y0^2))-y0*((x2-x1)-y0*(atan(x2/y0)-atan(x1/y0)));
    yConst11=-x1/2*y0*log((x2^2+y0^2)/(x1^2+y0^2))+y0*((x2-x1)-y0*(atan(x2/y0)-atan(x1/y0))); 
    B1=(4*NU-3)/2*(x2*(x2*log(x2^2+y0^2)-x1*log(x1^2+y0^2)...
        -2*((x2-x1)-y0*(atan(x2/y0)-atan(x1/y0))))-...
        1/2*((x2^2+y0^2)*log(x2^2+y0^2)-(x1^2+y0^2)*log(x1^2+y0^2)-(x2^2-x1^2)));
    B2=a^2*Treat2+b^2*y2Const2+2*a*b*yConst6;
    B3=c^2*Treat2+a^2*y2Const2+2*a*c*yConst6;
    B4=a*c*Treat2+a*b*y2Const2+(a^2+b*c)*yConst6;                    
    B5=(4*NU-3)/2*(-x1*(x2*log(x2^2+y0^2)-x1*log(x1^2+y0^2)...
        -2*((x2-x1)-y0*(atan(x2/y0)-atan(x1/y0))))+...
        1/2*((x2^2+y0^2)*log(x2^2+y0^2)-(x1^2+y0^2)*log(x1^2+y0^2)-(x2^2-x1^2)));
    B6=a^2*Treat6+b^2*y2Const7+2*a*b*yConst11;
    B7=c^2*Treat6+a^2*y2Const7+2*a*c*yConst11;
    B8=a*c*Treat6+a*b*y2Const7+(a^2+b*c)*yConst11;
    u11b=B1+B2+temp1;
    u22b=B1+B3+temp2;
    u12b=B4;
    u21b=u12b;
    u11a=B5+B6+temp1;
    u22a=B5+B7+temp2;
    u12a=B8;
    u21a=u12a;
    Gi=[u11b u12b u11a u12a;u21b u22b u21a u22a];
    Gi=Constp*double(Gi);
end   
%  Ui的求解
if  En1==np
    ersa=1-log(LEE);
    Hi=(1-2*NU)/4/pi/(1-NU)*[0 ersa 0 -1;-ersa 0 1 0];
    
elseif En2==np
    ersa=-1+log(LEE);
    Hi=(1-2*NU)/4/pi/(1-NU)*[0 1 0 ersa ;-1 0 -ersa 0 ];
else
    Constq=-LEE/4/pi/(1-NU)/(x2-x1)^2;
    %                     (7-8*NU)/2
    yConst2=x2*(atan(x2/y0)-atan(x1/y0))-y0*(log(x2^2+y0^2)-log(x1^2+y0^2))/2;
    yConst7=-x1*(atan(x2/y0)-atan(x1/y0))+y0*(log(x2^2+y0^2)-log(x1^2+y0^2))/2;
    yConst3=x2/2*((atan(x2/y0)-atan(x1/y0))-y0*(x2/(x2^2+y0^2)-x1/(x1^2+y0^2)))...
            -(y0*log((x2^2+y0^2)/(x1^2+y0^2))+y0^3*(1/(x2^2+y0^2)-1/(x1^2+y0^2)))/2;
    yConst8=-x1/2*((atan(x2/y0)-atan(x1/y0))-y0*(x2/(x2^2+y0^2)-x1/(x1^2+y0^2)))...
                 +(y0*log((x2^2+y0^2)/(x1^2+y0^2))+y0^3*(1/(x2^2+y0^2)-1/(x1^2+y0^2)))/2;
    y2Const5=(y0^2*(x2-x1)/(x1^2+y0^2)-y0*(atan(x2/y0)-atan(x1/y0)))/2;
    y2Const10=(y0^2*(x1-x2)/(x2^2+y0^2)+y0*(atan(x2/y0)-atan(x1/y0)))/2;
    y3Const4=x2/2*(y0*(x2/(x2^2+y0^2)-x1/(x1^2+y0^2))+(atan(x2/y0)-atan(x1/y0)))...
                  +y0^3*(1/(x2^2+y0^2)-1/(x1^2+y0^2))/2;  
    y3Const9=-x1/2*(y0*(x2/(x2^2+y0^2)-x1/(x1^2+y0^2))+(atan(x2/y0)-atan(x1/y0)))...
                    -y0^3*(1/(x2^2+y0^2)-1/(x1^2+y0^2))/2;
    Const6=x2/2*log((x2^2+y0^2)/(x1^2+y0^2))-((x2-x1)-y0*(atan(x2/y0)-atan(x1/y0))); 
    Const11=-x1/2*log((x2^2+y0^2)/(x1^2+y0^2))+(x2-x1)-y0*(atan(x2/y0)-atan(x1/y0)); 
    A1=(1-2*NU)*yConst2;
    A2=2*(a^2*yConst3+b^2*y3Const4+2*a*b*y2Const5);
    A3=2*(c^2*yConst3+a^2*y3Const4+2*a*c*y2Const5);
    A4=2*(a*c*yConst3+a*b*y3Const4+(a^2+b*c)*y2Const5);
    A5=-(1-2*NU)*((a*n2-c*n1)*Const6+(b*n2-a*n1)*yConst2);
    A6=(1-2*NU)*yConst7;
    A7=2*(a^2*yConst8+b^2*y3Const9+2*a*b*y2Const10);
    A8=2*(c^2*yConst8+a^2*y3Const9+2*a*c*y2Const10);
    A9=2*(a*c*yConst8+a*b*y3Const9+(a^2+b*c)*y2Const10);
    A10=-(1-2*NU)*((a*n2-c*n1)*Const11+(b*n2-a*n1)*yConst7);
    h11b=A1+A2;
    h22b=A1+A3;
    h12b=A4+A5;
    h21b=A4-A5;
    h11a=A6+A7;
    h22a=A6+A8;
    h12a=A9+A10;
    h21a=A9-A10;
    Hi=[h11b h12b h11a h12a;h21b h22b h21a h22a];
    Hi=Constq*double(Hi);
end
%%%计算体积力Bi
GravC=LEE/(x2-x1)/8/pi/SG;
if   En1==np
     GravC1=x2^2*log(x2^2)/2;
     GravC2=0;
elseif  En2==np
    GravC1=-x1^2*log(x1^2)/2;
    GravC2=0;
else
    GravC1=((x2^2+y0^2)*log(x2^2+y0^2)-(x1^2+y0^2)*log(x1^2+y0^2))/2;
    GravC2=y0*(x2*log(x2^2+y0^2)-x1*log(x1^2+y0^2)-2*((x2-x1)-y0*(atan(x2/y0)-atan(x1/y0)))+(x2-x1));
end   
    VForce=Den*AG*(n1*(a*GravC1+b*GravC2)+n2*(c*GravC1+a*GravC2))...
          -1/2/(1-NU)*Den*AG'*[(a*GravC1+b*GravC2) (c*GravC1+a*GravC2)]'*NE';
	VFi=VFi-VForce*GravC;

