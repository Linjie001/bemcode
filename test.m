        if ON==1
            ElementAiNode=[ElementAi(:,2) Node(ElementAi(:,3),2:3) Node(ElementAi(:,4),2:3)];
            iE=ElementAi(find(((ElementAiNode(:,2)-xp(1)).*(ElementAiNode(:,4)-xp(1))<=0).*...
            ((ElementAiNode(:,3)-xp(2)).*(ElementAiNode(:,5)-xp(2))<=0)),2);
            iE=iE(1);
            Nodeb=Element(iE,3);
            Nodea=Element(iE,4); 
            Lxpa=norm(xp'-Node(Nodea,2:3));
            Rb=Lxpa/LE(Element(iE,2));
            Ra=1-Rb;
            iNterUi=Ra*UPknow(Nodea,2:3)+Rb*UPknow(Nodeb,2:3);
            iNterU=cat(1,iNterU,[iA InNode(iKP,:) iNterUi]);
            
            PA=Ra*UPknow(Nodea,4:5)+Rb*UPknow(Nodeb,6:7);
            S11bar=PA*NE(iE,:)';
            S12bar=PA*DL(iE,:)';
            Strain22bar=diff(UPknow([Nodeb Nodea],2:3))*DL(iE,:)'/LE(iE);
            S22bar=1/(1-NU)*(NU*S11bar+2*SG*Strain22bar);
            
            MC=[NE(iE,1) NE(iE,2);-NE(iE,2) NE(iE,1)];
            Sxp=MC'*[S11bar S12bar;S12bar S22bar]*MC;

            iNterP=cat(1,iNterP,[iA InNode(iKP,:) [Sxp(1) Sxp(2) Sxp(4)]]);
            continue
       end  
       if IN==1
            AUi=zeros(2,1);
            APi=zeros(2,1);
            ADi=zeros(3,1);
            ASi=zeros(3,1);
            FVjj=zeros(2,1);
            SVjj=zeros(3,1);
        %   单元循环
            for iE=1:size(ElementAi,1)
                AD=zeros(3,4);
                AS=zeros(3,4);
                Gi=zeros(2,4);
                Hi=zeros(2,4);
                SV11=0;
                SV12=0;
                SV22=0;
                FV=0;
                ElNo=ElementAi(iE,2);
                En1=Element(ElNo,3);
                En2=Element(ElNo,4);
                %  需要进行局部坐标转换 Localxp=[0 0];
%                 定义的常数
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
                LocalNode(En1,:)=[1 (Trans*(Node(En1,2:3)'-xp))'];
                LocalNode(En2,:)=[2 (Trans*(Node(En2,2:3)'-xp))'];
                
                x2=LocalNode(En2,2);
                x1=LocalNode(En1,2);
                y0=LocalNode(En2,3);
% Ui
                Constp=LE(ElNo)/8/pi/SG/(1-NU)/(x2-x1)^2;
%                     (7-8*NU)/2
                temp1=-(7-8*NU)/2*(x2*(x2-x1)-(x2^2-x1^2)/2);
                temp2=-(7-8*NU)/2*(-x1*(x2-x1)+(x2^2-x1^2)/2);
%                     temp1=0;
%                     temp2=0
                Treat2=x2*((x2-x1)-y0*(atan(x2/y0)-atan(x1/y0)))...
                    -((x2^2-x1^2)/2-y0^2/2*log((x2^2+y0^2)/(x1^2+y0^2)));
                Treat6=-x1*((x2-x1)-y0*(atan(x2/y0)-atan(x1/y0)))...
                    +((x2^2-x1^2)/2-y0^2/2*log((x2^2+y0^2)/(x1^2+y0^2)));
                h1=x2*atan(x2/y0)-x1*atan(x1/y0)-y0/2*(log(x2^2+y0^2)-log(x1^2+y0^2));
                Const2=1/y0*(-(x2-x1)*atan(x1/y0)+h1);
                Const3=x2/2*(1/y0*(atan(x2/y0)-atan(x1/y0))-(x2/(x2^2+y0^2)-x1/(x1^2+y0^2)))...
                   -(log((x2^2+y0^2)/(x1^2+y0^2))+y0^2*(1/(x2^2+y0^2)-1/(x1^2+y0^2)))/2;
                Const4=x2/2/y0^2*((x2/(x2^2+y0^2)-x1/(x1^2+y0^2))+(atan(x2/y0)-atan(x1/y0))/y0)...
                      +(1/(x2^2+y0^2)-1/(x1^2+y0^2))/2;  
                Const5=((x2-x1)/(x1^2+y0^2)-(atan(x2/y0)-atan(x1/y0))/y0)/2;
                Const6=x2/2*log((x2^2+y0^2)/(x1^2+y0^2))...
                      -((x2-x1)-y0*(atan(x2/y0)-atan(x1/y0)));
                Const7=1/y0*((x2-x1)*atan(x2/y0)-h1);
                Const8=-x1/2*(1/y0*(atan(x2/y0)-atan(x1/y0))-(x2/(x2^2+y0^2)-x1/(x1^2+y0^2)))...
                       +(log((x2^2+y0^2)/(x1^2+y0^2))+y0^2*(1/(x2^2+y0^2)-1/(x1^2+y0^2)))/2;
                Const9=-x1/2/y0^2*((x2/(x2^2+y0^2)-x1/(x1^2+y0^2))+(atan(x2/y0)-atan(x1/y0))/y0)...
                      -(1/(x2^2+y0^2)-1/(x1^2+y0^2))/2;    
                Const10=((x1-x2)/(x2^2+y0^2)+(atan(x2/y0)-atan(x1/y0))/y0)/2;
                Const11=-x1/2*log((x2^2+y0^2)/(x1^2+y0^2))...
                   +(x2-x1)-y0*(atan(x2/y0)-atan(x1/y0)); 

                B1=(4*NU-3)/2*(x2*(x2*log(x2^2+y0^2)-x1*log(x1^2+y0^2)...
                    -2*((x2-x1)-y0*(atan(x2/y0)-atan(x1/y0))))-...
                    1/2*((x2^2+y0^2)*log(x2^2+y0^2)-(x1^2+y0^2)*log(x1^2+y0^2)-(x2^2-x1^2)));
                B2=a^2*Treat2+b^2*y0^2*Const2+2*a*b*y0*Const6;
                B3=c^2*Treat2+a^2*y0^2*Const2+2*a*c*y0*Const6;
                B4=a*c*Treat2+(a^2+b*c)*y0*Const6+a*b*y0^2*Const2;                    
                B5=(4*NU-3)/2*(-x1*(x2*log(x2^2+y0^2)-x1*log(x1^2+y0^2)...
                    -2*((x2-x1)-y0*(atan(x2/y0)-atan(x1/y0))))+...
                    1/2*((x2^2+y0^2)*log(x2^2+y0^2)-(x1^2+y0^2)*log(x1^2+y0^2)-(x2^2-x1^2)));
                B6=a^2*Treat6+b^2*y0^2*Const7+2*a*b*y0*Const11;
                B7=c^2*Treat6+a^2*y0^2*Const7+2*a*c*y0*Const11;
                B8=a*c*Treat6+(a^2+b*c)*y0*Const11+a*b*y0^2*Const7;

                u11b=B1+B2+temp1;
                u22b=B1+B3+temp2;
                u12b=B4;
                u21b=u12b;

                u11a=B5+B6+temp1;
                u22a=B5+B7+temp2;
                u12a=B8;
                u21a=u12a;
                Gi=[u11b u12b u11a u12a;u21b u22b u21a u22a];
                Gi=double(Gi);
% Hi;
                 Constq=-LE(ElNo)/4/pi/(1-NU)/(x2-x1)^2;
                 A1=(1-2*NU)*y0*Const2;
                 A2=2*y0*(a^2*Const3+b^2*y0^2*Const4+2*a*b*y0*Const5);
                 A3=2*y0*(c^2*Const3+a^2*y0^2*Const4+2*a*c*y0*Const5);
                 A4=2*y0*(a*c*Const3+a*b*y0^2*Const4+(a^2+b*c)*y0*Const5);
                 A5=-(1-2*NU)*((a*NE(ElNo,2)-c*NE(ElNo,1))*Const6+(b*NE(ElNo,2)-a*NE(ElNo,1))*y0*Const2);

                 h11b=A1+A2;
                 h22b=A1+A3;
                 h12b=A4+A5;
                 h21b=A4-A5;
                 A6=(1-2*NU)*y0*Const7;
                 A7=2*y0*(a^2*Const8+b^2*y0^2*Const9+2*a*b*y0*Const10);
                 A8=2*y0*(c^2*Const8+a^2*y0^2*Const9+2*a*c*y0*Const10);
                 A9=2*y0*(a*c*Const8+a*b*y0^2*Const9+(a^2+b*c)*y0*Const10);
                 A10=-(1-2*NU)*((a*NE(ElNo,2)-c*NE(ElNo,1))*Const11+(b*NE(ElNo,2)-a*NE(ElNo,1))*y0*Const7);

                 h11a=A6+A7;
                 h22a=A6+A8;
                 h12a=A9+A10;
                 h21a=A9-A10;
                 Hi=[h11b h12b h11a h12a;h21b h22b h21a h22a];
                 Hi=double(Hi);
%                  --求物体内部点的应力所需要的变量
                Const=LE(ElNo)/4/pi/(1-NU)/(x2-x1)^2;
                rou1=a*Const6+b*y0*Const2;
                je1=Const6-y0^2*Const5;
                rou2=a^3*je1+3*a^2*b*y0*Const3+3*a*b^2*y0^2*Const5+b^3*y0^3*Const4;
                rou3=a*Const11+b*y0*Const7;
                je2=Const11-y0^2*Const10;
                rou4=a^3*je2+3*a^2*b*y0*Const8+3*a*b^2*y0^2*Const10+b^3*y0^3*Const9;
                rou5=-c*Const6-a*y0*Const2;
                rou6=a^2*c*je1+(a^3+2*a*b*c)*y0*Const3+(b^2*c+2*a^2*b)*y0^2*Const5+a*b^2*y0^3*Const4;
                rou7=-c*Const11-a*y0*Const7;
                rou8=a^2*c*je2+(a^3+2*a*b*c)*y0*Const8+(b^2*c+2*a^2*b)*y0^2*Const10+a*b^2*y0^3*Const9;
                rou9=a*c^2*je1+(2*a^2*c+b*c^2)*y0*Const3+(a^3+2*a*b*c)*y0^2*Const5+a^2*b*y0^3*Const4;
                rou10=a*c^2*je2+(2*a^2*c+b*c^2)*y0*Const8+(a^3+2*a*b*c)*y0^2*Const10+a^2*b*y0^3*Const9;
                rou11=c^3*je1+3*c^2*a*y0*Const3+3*c*a^2*y0^2*Const5+a^3*y0^3*Const4;
                rou12=c^3*je2+3*c^2*a*y0*Const8+3*c*a^2*y0^2*Const10+a^3*y0^3*Const9;
                D111q1=(1-2*NU)*rou1+2*rou2;
                D111q2=(1-2*NU)*rou3+2*rou4; 
                D211q1=(1-2*NU)*rou5+2*rou6;
                D211q2=(1-2*NU)*rou7+2*rou8;
                D112q1=-(1-2*NU)*rou5+2*rou6;
                D112q2=-(1-2*NU)*rou7+2*rou8;
                D212q1=(1-2*NU)*rou1+2*rou9;
                D212q2=(1-2*NU)*rou3+2*rou10;
                D122q1=-(1-2*NU)*rou1+2*rou9;
                D122q2=-(1-2*NU)*rou3+2*rou10;
                D222q1=-(1-2*NU)*rou5+2*rou11;
                D222q2=-(1-2*NU)*rou7+2*rou12;
                AD=Const*[D111q1 D211q1 D111q2 D211q2;D112q1 D212q1 D112q2 D212q2;D122q1 D222q1 D122q2 D222q2];
%                 ---ASi-
               Const=LE(ElNo)*SG/2/pi/(1-NU)/(x2-x1)^2;
               h2=1/4/y0^2*(x2/(x2^2+y0^2)^2-x1/(x1^2+y0^2)^2+...
                  3/2/y0^2*(x2/(x2^2+y0^2)-x1/(x1^2+y0^2)+1/y0*(atan(x2/y0)-atan(x1/y0))));
               cst1=1/2*(1/(x1^2+y0^2)-1/(x2^2+y0^2)+y0^2/2*(1/(x2^2+y0^2)^2-1/(x1^2+y0^2)^2));
               cst2=1/2/y0^2*(x2/(x2^2+y0^2)-x1/(x1^2+y0^2)+1/y0*(atan(x2/y0)-atan(x1/y0)))-y0^2*h2;
               cst3=-1/4*(1/(x2^2+y0^2)^2-1/(x1^2+y0^2)^2);
               cst4=1/y0*(atan(x2/y0)-atan(x1/y0))-2*y0^2*cst2-y0^4*h2;
               wer1=a*Const5+b*y0*Const4;
               wer2=a*Const10+b*y0*Const9;
               wer3=c*Const5+a*y0*Const4;
               wer4=c*Const10+a*y0*Const9;
               mj1=a^3*(x2*cst1-cst4)+3*a^2*b*y0*(x2*cst2-cst1)+3*a*b^2*y0^2*(x2*cst3-cst2)+b^3*y0^3*(x2*h2-cst3);
               mj2=a^3*(cst4-x1*cst1)+3*a^2*b*y0*(cst1-x1*cst2)+3*a*b^2*y0^2*(cst2-x1*cst3)+b^3*y0^3*(cst3-x1*h2);
               mj3=a^2*c*(x2*cst1-cst4)+(a^3+2*a*b*c)*y0*(x2*cst2-cst1)+(b^2*c+2*a^2*b)*y0^2*(x2*cst3-cst2)+b^2*a*y0^3*(x2*h2-cst3);
               mj4=a^2*c*(cst4-x1*cst1)+(a^3+2*a*b*c)*y0*(cst1-x1*cst2)+(b^2*c+2*a^2*b)*y0^2*(cst2-x1*cst3)+b^2*a*y0^3*(cst3-x1*h2);
               mj5=a*c^2*(x2*cst1-cst4)+(2*a^2*c+b*c^2)*y0*(x2*cst2-cst1)+(a^3+2*a*b*c)*y0^2*(x2*cst3-cst2)+a^2*b*y0^3*(x2*h2-cst3);
               mj6=a*c^2*(cst4-x1*cst1)+(2*a^2*c+b*c^2)*y0*(cst1-x1*cst2)+(a^3+2*a*b*c)*y0^2*(cst2-x1*cst3)+a^2*b*y0^3*(cst3-x1*h2);
               mj7=c^3*(x2*cst1-cst4)+3*c^2*a*y0*(x2*cst2-cst1)+3*c*a^2*y0^2*(x2*cst3-cst2)+a^3*y0^3*(x2*h2-cst3);
               mj8=c^3*(cst4-x1*cst1)+3*c^2*a*y0*(cst1-x1*cst2)+3*c*a^2*y0^2*(cst2-x1*cst3)+a^3*y0^3*(cst3-x1*h2);
               me1=wer1-4*mj1;
               me2=a^2*Const3+2*a*b*y0*Const5+b^2*y0^2*Const4;
               me3=wer2-4*mj2;
               me4=a^2*Const8+2*a*b*y0*Const10+b^2*y0^2*Const9;
               me5=(1-2*NU)*wer3-4*mj3;
               me6=a*c*Const3+(a^2+b*c)*y0*Const5+a*b*y0^2*Const4;
               me7=(1-2*NU)*wer4-4*mj4;
               me8=a*c*Const8+(a^2+b*c)*y0*Const10+a*b*y0^2*Const9;
               me9=NU*wer3-4*mj3;
               me10=NU*wer4-4*mj4;
               me11=NU*wer1-4*mj5;
               me12=c^2*Const3+2*a*c*y0*Const5+a^2*y0^2*Const4;
               me13=NU*wer2-4*mj6;
               me14=c^2*Const8+2*a*c*y0*Const10+a^2*y0^2*Const9;
               me15=(1-2*NU)*wer1-4*mj5;
               me16=(1-2*NU)*wer2-4*mj6;
               me17=wer3-4*mj7;
               me18=wer4-4*mj8;
               n1=NE(ElNo,1);
               n2=NE(ElNo,2);
               S111q1=2*y0*me1+2*n1*me2+n1*Const2;
               S111q2=2*y0*me3+2*n1*me4+n1*Const7;
               S211q1=2*y0*me5+4*NU*n1*me6+(1-2*NU)*2*n2*me2-(1-4*NU)*n2*Const2;
               S211q2=2*y0*me7+4*NU*n1*me8+(1-2*NU)*2*n2*me4-(1-4*NU)*n2*Const7;
               S112q1=2*y0*me9+(2-2*NU)*n1*me6+2*NU*n2*me2+(1-2*NU)*n2*Const2;
               S112q2=2*y0*me10+(2-2*NU)*n1*me8+2*NU*n2*me4+(1-2*NU)*n2*Const7;
               S212q1=2*y0*me11+2*n1*NU*me12+(2-2*NU)*n2*me6+(1-2*NU)*n1*Const2;
               S212q2=2*y0*me13+2*n1*NU*me14+(2-2*NU)*n2*me8+(1-2*NU)*n1*Const7;
               S122q1=2*y0*me15+4*NU*n2*me6+(1-2*NU)*2*n1*me12-(1-4*NU)*n1*Const2;
               S122q2=2*y0*me16+4*NU*n2*me8+(1-2*NU)*2*n1*me14-(1-4*NU)*n1*Const7;
               S222q1=2*y0*me17+2*n2*me12+n2*Const2;
               S222q2=2*y0*me18+2*n2*me14+n2*Const7;
               AS=Const*[S111q1 S211q1 S111q2 S211q2;S112q1 S212q1 S112q2 S212q2;S122q1 S222q1 S122q2 S222q2];      
                           
%                 ----------------
               AUi(:,1)=AUi(:,1)+Gi*Constp*[UPknow(En1,6:7) UPknow(En2,4:5)]';
               APi(:,1)=APi(:,1)+Hi*Constq*[UPknow(En1,2:3) UPknow(En2,2:3)]';
               ADi(:,1)= ADi(:,1)+AD*[UPknow(En1,6:7) UPknow(En2,4:5)]';
               ASi(:,1)= ASi(:,1)+AS*[UPknow(En1,2:3) UPknow(En2,2:3)]';
               FVjj(:,1)=FVjj(:,1)+FV;
               SVjj(:,1)=SVjj(:,1)+[SV11 SV12 SV22]'; 
            end 