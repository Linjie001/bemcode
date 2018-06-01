function [H G VF]=HGFA(np,ElementAi,IntMethod)
global Element GaussXW Node LE MatiA Mat AG ElNoC
% PQNM=0采用解析解
% PQNM=1采用数值解
PQNM=1;
RIndex=1:2;
H=zeros(2,size(Node,1)*2);
G=zeros(2,size(Node,1)*4);
VF=zeros(2,1);
    for iE=1:size(ElementAi,1)
             %求H矩阵
            ElNo=ElementAi(iE,2);
            En1=Element(ElNo,3);
            En2=Element(ElNo,4);
            
            if IntMethod==0
                xl=([1-x x+1]/2*Node([En1 En2],2:end)).';
                r=xl-xp;
                rl=sqrt(r.'*r);
                drd1=r(1)/rl;
                drd2=r(2)/rl;
                drdn=[drd1 drd2]*NE(ElNo,:)';
                NShape=[(1-x)/2 0 (1+x)/2 0;0 (1-x)/2 0 (1+x)/2];
                p11=drdn*((1-2*NU)+2*drd1^2);
                p12=drdn*(2*drd1*drd2)-(1-2*NU)*(drd1*NE(ElNo,2)-drd2*NE(ElNo,1));
                p21=drdn*(2*drd1*drd2)+(1-2*NU)*(drd1*NE(ElNo,2)-drd2*NE(ElNo,1));
                p22=drdn*((1-2*NU)+2*drd2^2); 
                P=-1/4/pi/(1-NU)/rl*[p11 p12;p21 p22];
                FPN=P*NShape*LE(ElNo)/2;
                if (En1==np | En2==np)
                    if En1==np
                        ElNoC=ElementAi(find(ElementAi(:,4)==np),2);
                    else
                        ElNoC=ElementAi(find(ElementAi(:,3)==np),2);
                    end
                    if En1==np %单元首节点为源点
                        Hi(1:2,3:4)=quadv(inline(FPN(:,3:4)),-1,1);
                        Hi(1:2,1:2)=(1-2*NU)/4/pi/(1-NU)*log(LE(ElNoC)/LE(ElNo))*[0 1;-1 0];
                    else                    
                        Hi(1:2,3:4)=0;
                        Hi(1:2,1:2)=quadv(inline(FPN(:,1:2)),-1,1);
                    end 
                else
                    Hi=quadv(inline(FPN),-1,1);
                end           
                Hi=double(Hi);
                CIndex=[2*En1+[-1 0],2*En2+[-1 0]];
                H(RIndex,CIndex)=H(RIndex,CIndex)+Hi;
                %求G矩阵
                if (En1==np | En2==np)
                    AX=diff(Node([En1 En2],2))/2;
                    AY=diff(Node([En1 En2],3))/2;
                    SR=(AX^2+AY^2)^0.5;  
                    DE=4*pi*SG*(1-NU);               
                    SR2=2*SR;   %一个单元的长度
                    DE2=2*DE;
                    XXU=3-4*NU;
                    X21=2*AX;
                    Y21=2*AY;
                    AG(1)=(SR*XXU*(1.5-log(SR2))+X21^2/(2*SR2));
                    AG(4)=(SR*XXU*(1.5-log(SR2))+Y21^2/(2*SR2));
                    BG(1)=(SR*XXU*(0.5-log(SR2))+X21^2/(2*SR2));
                    BG(4)=(SR*XXU*(0.5-log(SR2))+Y21^2/(2*SR2));
                    AG(2)=X21*Y21/(2*SR2);
                    BG(2)=AG(2);
                    AG(3)=AG(2);
                    BG(3)=BG(2);
                    Gi = [AG(1) AG(2) BG(1) BG(2);AG(3) AG(4) BG(3) BG(4)]/DE2                    
                else
                    %need to be added                    
                end                
                CIndex=[4*En1+[-1 0] 4*En2+[-3 -2]];
                G(RIndex,CIndex)=G(RIndex,CIndex)+Gi;
                %%%计算体积力Bi
                VForce=rl*(2*log(rl)+1)*(Den*AG*([drd1 drd2]*NE(ElNo,:)')-...
                                         1/2/(1-NU)*Den*AG'*[drd1 drd2]'*NE(ElNo,:)');
                VForce=quadv(inline(VForce),-1,1)*LE(ElNo)/2;
                VForce=-VForce/8/pi/SG;
                VF(RIndex)=VF(RIndex)+VForce;
                
            elseif IntMethod==1 | IntMethod==2
                syms x
                xp=Node(np,2:end)'
                xl=([1-x x+1]/2*Node([En1 En2],2:end)).'; 
                r=xl-xp;
                rl=sqrt(r.'*r);            
                drd1=r(1)/rl;
                drd2=r(2)/rl;
                drdn=[drd1 drd2]*NE(ElNo,:)';
                NShape=[(1-x)/2 0 (1+x)/2 0;0 (1-x)/2 0 (1+x)/2];
                p11=drdn*((1-2*NU)+2*drd1^2);
                p12=drdn*(2*drd1*drd2)-(1-2*NU)*(drd1*NE(ElNo,2)-drd2*NE(ElNo,1));
                p21=drdn*(2*drd1*drd2)+(1-2*NU)*(drd1*NE(ElNo,2)-drd2*NE(ElNo,1));
                p22=drdn*((1-2*NU)+2*drd2^2); 
                P=-1/4/pi/(1-NU)/rl*[p11 p12;p21 p22];
                FPN=P*NShape*LE(ElNo)/2;
                if (En1==np | En2==np)
                    if En1==np
                        ElNoC=ElementAi(find(ElementAi(:,4)==np),2);
                    else
                        ElNoC=ElementAi(find(ElementAi(:,3)==np),2);
                    end
                    if En1==np %单元首节点为源点
                        if IntMethod==1
                            Hi(1:2,3:4)=int(FPN(:,3:4),-1,1);
                        else
                            Hi(1:2,3:4)=quadv(inline(FPN(:,3:4)),-1,1);
                        end                        
                        Hi(1:2,1:2)=(1-2*NU)/4/pi/(1-NU)*log(LE(ElNoC)/LE(ElNo))*[0 1;-1 0];
                    else                    
                        Hi(1:2,3:4)=0;
                        if IntMethod==1
                            Hi(1:2,1:2)=int(FPN(:,1:2),-1,1);
                        else
                            Hi(1:2,1:2)=quadv(inline(FPN(:,1:2)),-1,1);
                        end
                    end 
                else
                    if IntMethod==1
                        Hi=int(FPN,-1,1);
                    else
                        Hi=quadv(inline(FPN),-1,1);
                    end
                end           
                Hi=double(Hi);
                CIndex=[2*En1+[-1 0],2*En2+[-1 0]];
                H(RIndex,CIndex)=H(RIndex,CIndex)+Hi;
                %求G矩阵
                u11=(4*NU-3)*log(rl)+drd1^2-(7-8*NU)/2;
                u12=drd1*drd2;
                u21=u12;
                u22=(4*NU-3)*log(rl)+drd2^2-(7-8*NU)/2;
                U=[u11 u12;u21 u22];
                Gi=U*NShape;
                if IntMethod==1
                    Gi=int(Ui,-1,1);
                else
                    Gi=quadv(inline(Gi),-1,1);
                end
                Gi=double(Gi);
                Gi=1/8/pi/SG/(1-NU)*Gi*LE(ElNo)/2;
                CIndex=[4*En1+[-1 0] 4*En2+[-3 -2]];
                G(RIndex,CIndex)=G(RIndex,CIndex)+Gi;
                %%%计算体积力Bi
                VForce=rl*(2*log(rl)+1)*(Den*AG*([drd1 drd2]*NE(ElNo,:)')-...
                                         1/2/(1-NU)*Den*AG'*[drd1 drd2]'*NE(ElNo,:)');
                if IntMethod==1
                    VForce=int(VForce,-1,1)*LE(ElNo)/2;
                else
                    VForce=quadv(inline(VForce),-1,1)*LE(ElNo)/2;
                end
                VForce=double(VForce);
                VForce=-VForce/8/pi/SG;
                VF(RIndex)=VF(RIndex)+VForce;                
            elseif IntMethod==3
                if En1==np
                    ElNoC=ElementAi(find(ElementAi(:,4)==np),2);
                elseif En2==np
                    ElNoC=ElementAi(find(ElementAi(:,3)==np),2);
                end
                if PQNM==1
                     [Hi Gi VFi]=HGFE(GaussXW,Node([np En1 En2],:),LE(ElNoC),Mat(MatiA,:),AG);
                else
%                     解析解
                    [Hi Gi VFi]=ModiHGFE(GaussXW,Node([np En1 En2],:),LE(ElNoC),Mat(MatiA,:),AG);
                end
                CIndex=[2*En1+[-1 0],2*En2+[-1 0]];
                H(RIndex,CIndex)=H(RIndex,CIndex)+Hi;
                CIndex=[4*En1+[-1 0] 4*En2+[-3 -2]];
                G(RIndex,CIndex)=G(RIndex,CIndex)+Gi;
                VF(RIndex)=VF(RIndex)+VFi;
            end
        end