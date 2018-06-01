function [Hi Gi VFi]=HGF(GaussXW,Node,LEC,Mat,AG)
% ENCo,xp,GaussXW,GaussPt
    np=Node(1,1);	En1=Node(2,1);  En2=Node(3,1);    xp=Node(1,2:3)';    
	Den=Mat(1);    NU=Mat(3);    SG=Mat(4);    
    Hi=zeros(2,4);	Gi=zeros(2,4);	VFi=zeros(2,1);
    
     DLE=Node(3,2:3)-Node(2,2:3);
     LEE=sum(DLE.*DLE,2).^0.5;
     DLE=DLE./[LEE LEE];
     NEE=DLE*[0 -1;1 0];
    
	for iG=1:size(GaussXW,2)
        x=GaussXW(1,iG);
        xl=([1-x x+1]/2*Node(2:3,2:3)).'; 
        r=xl-xp;
        rl=sqrt(r.'*r);            
        drd1=r(1)/rl;
        drd2=r(2)/rl;
        drdn=[drd1 drd2]*NEE';
        NShape=[(1-x)/2 0 (1+x)/2 0;0 (1-x)/2 0 (1+x)/2];
        %求H矩阵
        p11=drdn*((1-2*NU)+2*drd1^2);
        p12=drdn*(2*drd1*drd2)-(1-2*NU)*(drd1*NEE(2)-drd2*NEE(1));
        p21=drdn*(2*drd1*drd2)+(1-2*NU)*(drd1*NEE(2)-drd2*NEE(1));
        p22=drdn*((1-2*NU)+2*drd2^2); 
        P=-1/4/pi/(1-NU)/rl*[p11 p12;p21 p22];
        FPN=P*NShape*LEE/2;
        if (En1==np | En2==np)
            if En1==np %单元首节点为源点
                Hi(:,3:4)=Hi(:,3:4)+FPN(:,3:4)*GaussXW(2,iG);
                if iG==1
                    Hi(:,1:2)=(1-2*NU)/4/pi/(1-NU)*log(LEC/LEE)*[0 1;-1 0];
                end
            else                    
                Hi(:,3:4)=0;
                Hi(:,1:2)=Hi(:,1:2)+FPN(:,1:2)*GaussXW(2,iG);
            end 
        else
            Hi=Hi+FPN*GaussXW(2,iG);
        end
        %求G矩阵
        if En1~=np & En2~=np
            u11=(4*NU-3)*log(rl)+drd1^2-(7-8*NU)/2;
            u12=drd1*drd2;
            u21=u12;
            u22=(4*NU-3)*log(rl)+drd2^2-(7-8*NU)/2;
            U=1/8/pi/SG/(1-NU)*LEE/2*[u11 u12;u21 u22];
            Gi=Gi+U*NShape*GaussXW(2,iG); 
        elseif  iG==1 & (En1==np | En2==np)
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
	end
	 %%%计算体积力Bi
	VForce=rl*(2*log(rl)+1)*(Den*AG*([drd1 drd2]*NEE')-...
                            1/2/(1-NU)*Den*AG'*[drd1 drd2]'*NEE');
	VFi=VFi-VForce*GaussXW(2,iG)*LEE/2/8/pi/SG;                              
end

