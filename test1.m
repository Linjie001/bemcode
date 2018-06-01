% 该程序在求解内部点位移时对单元采用解析解
clc
clear
global Element GaussXW Node LE MatiA Mat AG  ElNoC
[ProblemType,Mat,Node,Area,Element,BCU,InNode]=Input('data.txt');
if ProblemType==2
    Mat(:,2)=Mat(:,2)*(1+2*Mat(:,3))/(1+Mat(:,3))^2;
    Mat(:,3)=Mat(:,3)/(1+Mat(:,3));    
end
AG=[0 -9.8]';       %   acceleration of gravity
IntMethod=3;      %0:解析解   1:int    2:quadv   3:高斯积分
GaussPt=5;
MChau=0
% 以下是程序部分
set(0,'FormatSpacing','compact')
%前处理中公共边界的处理
%公共边界的查找
CommonElement=zeros(0,4);
for iE=1:size(Element,1)
    IP=find(Element((iE+1):end,4)==Element(iE,3) &...
        Element((iE+1):end,3)==Element(iE,4));
    if ~isempty(IP)
       iE1=Element(iE,2);
       MatA1=Element(iE,1);
       MatA2=Element(iE+IP,1);
       if MatA1<MatA2
           CommonElement=cat(1,CommonElement,[iE1 MatA1 0 MatA2]);
       else
           CommonElement=cat(1,CommonElement,[iE1 MatA2 0 MatA1]);
       end
    end
end
%增加节点，修改对应单元中的节点号并增加其对应的边界条件

Element(:,end+2)=0;
for iA0=1:(size(Area,1)-1)
    for iA1=(iA0+1):size(Area,1)
        IndexComEl=find(CommonElement(:,2)==iA0 & CommonElement(:,4)==iA1);
        if isempty(IndexComEl)
            continue
        end
        for iCE=1:length(IndexComEl)
            E0=CommonElement(IndexComEl(iCE),1);
            E1=max(Element(:,2))+1;
            InE1=find(Element(:,2)==CommonElement(IndexComEl(iCE),1) & ...
            Element(:,1)==CommonElement(IndexComEl(iCE),4));
            Element(InE1,2)=E1;
            CommonElement(IndexComEl(iCE),3)=E1;
            InE0=find(Element(:,2)==E0);
            for iN = 1:2
                N0=Element(InE0,iN+2);
                N1=Element(InE1,5-iN);
%                李春光的程序
                inENa=find(Element(:,3)==N1 & Element(:,1)>iA0 & Element(:,13)==0);
                inENb=find(Element(:,4)==N1 & Element(:,1)>iA0 & Element(:,14)==0);
%                 邓琴的程序
%                 temp=iA0;
%                 if Element(InE0,13:14)==0 & Element(InE1,13:14)~=0
%                     temp=Element(InE1,1);
%                 end
%                 inENa=find(Element(:,3)==N1 & Element(:,1)~=temp & Element(:,13)==0);
%                 inENb=find(Element(:,4)==N1 & Element(:,1)~=temp & Element(:,14)==0);
                if  isempty(union(inENa,inENb)) & N0~=N1
                    continue
                end
                N11=size(Node,1)+1;
                
                if isempty(union(inENa,inENb)) %前面大面号，后面小面号
                    Node(N11,:)=[N11 Node(N0,2:3)];
                    inENa=find(Element(:,3)==N0 & Element(:,1)<iA1 & Element(:,13)==0);
                    inENb=find(Element(:,4)==N0 & Element(:,1)<iA1 & Element(:,14)==0);
                    Element(inENa,3)=N11;
                    Element(inENb,4)=N11;
                    Element(find(Element(:,3)==N11 & Element(:,1)==iA0),13)=1;
                    Element(find(Element(:,4)==N11 & Element(:,1)==iA0),14)=1;
                else    %前面小面号，后面大面号
                    Node(N11,:)=[N11 Node(N1,2:3)];
                    Element(inENa,3)=N11;
                    Element(inENb,4)=N11;
                    Element(find(Element(:,3)==N11 & Element(:,1)==iA1),13)=1;
                    Element(find(Element(:,4)==N11 & Element(:,1)==iA1),14)=1;
                end                
            end
            Element(InE0,5:8)=0;
            Element(InE1,5:8)=0;
        end        
    end
end
Element(:,13:end)=[];
%对单元进行排序
 [B,IX]=sort(Element(:,2));
 Element=Element(IX,:);
 
CommonNodes=zeros(0,2);  %[主节点 从节点]
for iE=1:size(CommonElement,1)
    for iN=1:2
        Na=Element(CommonElement(iE,1),iN+2);
        Nb=Element(CommonElement(iE,3),5-iN);
        if isempty(find(CommonNodes(:,1)==Na & CommonNodes(:,2)==Nb))
            if isempty(find(CommonNodes(:,1)==Na))  %每个节点在主节点只出现一次
                CommonNodes=cat(1,CommonNodes,[Na Nb]);
            else
                CommonNodes=cat(1,CommonNodes,[Nb Na]);
            end
        end
    end
end
%去除多余的公共点对,指该点四周全被各个面包围的点，不去除会增加多余的约束
RowTobeDel=zeros(0,1);
CommonNodes(:,end+1)=0;
for iE=1:size(CommonNodes,1)
    UsedRow=[];
    if CommonNodes(iE,end)==1
        continue
    end
    Node0=CommonNodes(iE,1);
    Nodei=CommonNodes(iE,2);
    Nodej=0;
    Row = find(CommonNodes(:,1)==Nodei);
    if CommonNodes(Row,2)~=Node0
        Nodej=CommonNodes(Row,2);
        UsedRow=cat(1,UsedRow,Row);
    end
    Row = find(CommonNodes(:,2)==Nodei);
    if CommonNodes(Row,1)~=Node0
        Nodej=CommonNodes(Row,1);
        UsedRow=cat(1,UsedRow,Row);
    end
    while  Nodej~=0        
        Row=find((CommonNodes(:,1)==Nodej) & (CommonNodes(:,2)~=Nodei)...
                |(CommonNodes(:,2)==Nodej) & (CommonNodes(:,1)~=Nodei));
        if isempty(Row)
            break
        end
        CommonNodes(Row,end)=1;
        Nodei=Nodej;
        Nodej=sum(CommonNodes(Row,1:2))-Nodei;
        if Nodej==Node0
            RowTobeDel=cat(1,RowTobeDel,Row);
            break
        end            
    end 
    CommonNodes(iE,end)=1;
    CommonNodes(UsedRow,end)=1;
end
CommonNodes(RowTobeDel,:)=[];
CommonNodes(:,end)=[];
% 求单元的法向
 DL=Node(Element(:,4),2:end)-Node(Element(:,3),2:end);
 LE=sum(DL.*DL,2).^0.5;
 DL=DL./[LE LE];
 NE=DL*[0 -1;1 0];
%边界节点已知位移量个数
 UPknow=zeros(size(Node,1),7);
 UPknow(:,1)=1:size(Node,1);
 IndexUPknow=zeros(size(Node,1),8);
 IndexUPknow(:,1)=1:size(Node,1);
 IUX=find(BCU(:,2)==1);
 IUY=find(BCU(:,2)==2);
%  IUXY=find(BCU(:,2)==3);
 
 UPknow(BCU(IUX),2)=BCU(IUX,3);
 UPknow(BCU(IUY),3)=BCU(IUY,3);
%  UPknow(BCU(IUXY),2:3)=BCU(IUXY,3)*[1 1]

 IndexUPknow(BCU(IUX),2)=1;
 IndexUPknow(BCU(IUY),3)=1;
%  IndexUPknow(BCU(IUXY),2:3)=ones(size(IUXY,1),2)
 
 
IPXb=Element(find(Element(:,7)==1),4);
IPYb=Element(find(Element(:,8)==1),4);
% IPXYb=find(BCP(:,2)==3)
IPXa=Element(find(Element(:,5)==1),3);
IPYa=Element(find(Element(:,6)==1),3);
% IPXYa=find(BCP(:,2)==13)

UPknow(IPXb,4)=Element(find(Element(:,7)==1),11);
UPknow(IPYb,5)=Element(find(Element(:,8)==1),12);
% UPknow(BCP(IPXYb),4:5)=BCP(IPXYb,3)*[1 1]
UPknow(IPXa,6)=Element(find(Element(:,5)==1),9);
UPknow(IPYa,7)=Element(find(Element(:,6)==1),10);
% UPknow(BCP(IPXYa),6:7)=BCP(IPXYa,3)*[1 1]

IndexUPknow(IPXb,4)=1;
IndexUPknow(IPYb,5)=1;
% IndexUPknow(BCP(IPXYb),4:5)=ones(size(BCP(IPXYb,3)))*[1 1]
IndexUPknow(IPXa,6)=1;
IndexUPknow(IPYa,7)=1;
% IndexUPknow(BCP(IPXYa),6:7)=ones(size(BCP(IPXYa,3)))*[1 1]

IndexUPknow(:,end)=sum(IndexUPknow(:,2:7),2);
IndexUPknow(:,end+1)=IndexUPknow(:,end);

%IndexU1已知位移边界，IndexU0未知位移边界

U1=cat(1,[BCU(IUX,1)*2-1,BCU(IUX,3)],[BCU(IUY,1)*2,BCU(IUY,3)]);
[IndexU1,IX] = sort(U1(:,1));
U1=U1(IX,2);
IndexU0=setdiff(1:size(Node,1)*2, IndexU1)';

%IndexP1已知应力边界，IndexP0未知应力边界
 P1=cat(1,[IPXb*4-3,UPknow(IPXb,4)],[IPYb*4-2,UPknow(IPYb,5)],...
[IPXa*4-1,UPknow(IPXa,6)],[IPYa*4,UPknow(IPYa,7)]);
 [IndexP1,IX] = sort(P1(:,1));
 P1=P1(IX,2);
IndexP0=setdiff(1:size(Node,1)*4, IndexP1)';

RowNos=size(Element,1)*2;
H=zeros(RowNos);
C=H;
G=zeros(RowNos,RowNos*2);
VF=zeros(RowNos,1);

if IntMethod==0
    syms x real
elseif IntMethod==1 | IntMethod==2
    syms x real
elseif IntMethod==3
    syms x;         %计算高斯点坐标及系数
    maple('with','orthopoly');
    f=maple('P',GaussPt,x);		%n为阶数
    c=sym2poly(f);
    x = sort(roots(c));
    p = legendre(GaussPt-1,x);
    c = 2*(1-x.^2)./(GaussPt*p(1,:)').^2;
    GaussXW = [x c]';
end
tic
for iA=1:max(Element(:,1))
    MatiA=Area(iA,2);
    Den=Mat(MatiA,1);
    NU=Mat(MatiA,3);
    SG=Mat(MatiA,4);
    ElementAi=Element(find(Element(:,1)==iA),:);
    NodeAi=ElementAi(:,3:4);
    NodeAi=reshape(NodeAi,[],1);
    NodeAi=union(NodeAi,NodeAi);
    
    for iN=1:size(NodeAi,1)
        disp(sprintf('Area=%d        iN=%d',iA,NodeAi(iN)))
        %求C矩阵 
        np=NodeAi(iN); 
        n1=ElementAi(find(ElementAi(:,4)==np),3);           
        n2=ElementAi(find(ElementAi(:,3)==np),4);
        a1=DirectionAngle(Node(np,2:end),Node(n1,2:end));
        a2=DirectionAngle(Node(np,2:end),Node(n2,2:end));
        if a2<a1
            a2=a2+2*pi;
        end
        da=a2-a1;
        aa=(a1+a2);
        nu2=2*(1-NU);
        Ci=eye(2)-1/2/pi/nu2*[nu2*da+cos(aa)*sin(da) sin(aa)*sin(da);...
                              sin(aa)*sin(da) nu2*da-cos(aa)*sin(da)];
        RIndex=2*np+[-1 0];
        C(RIndex,RIndex)=Ci;
        
        xp=Node(np,2:end)';
        [Hnp Gnp VFnp]=HGFA(np,ElementAi,IntMethod);
   
        H(RIndex,:)=H(RIndex,:)+Hnp;
        G(RIndex,:)=G(RIndex,:)+Gnp;
        VF(RIndex)=VF(RIndex)+VFnp;
    end
end

Gm=mean(Mat(:,4));
HH=C+H;
G=G*Gm;
P1=P1/Gm;
%计算
HGL=cat(2,HH(:,IndexU0),-G(:,IndexP0));
HGR=cat(2,-HH(:,IndexU1),G(:,IndexP1));
UP1=cat(1,U1,P1);
F=HGR*UP1;
F=F+VF;

%施加位移连续条件
for iN=1:size(CommonNodes,1)
    iN0=CommonNodes(iN,1);    %主节点号
    iN1=CommonNodes(iN,2);    %从节点号
    for iDOF=1:2             
        ColPos1=numel(IndexUPknow(1:iN1-1,2:3))-sum(sum(IndexUPknow(1:iN1-1,2:3)))+...
            1+(1-IndexUPknow(iN1,2))*(iDOF-1);           
        if IndexUPknow(iN0,iDOF+1)==0        %位移未知
            ColPos0=numel(IndexUPknow(1:iN0-1,2:3))-sum(sum(IndexUPknow(1:iN0-1,2:3)))+...
            	1+(1-IndexUPknow(iN0,2))*(iDOF-1);
        	HGL(end+1,[ColPos0 ColPos1])=[1 -1];
            F(end+1)=0;
        else
            HGL(end+1,ColPos1)=1;
            F(end+1)=BCU(find(BCU(:,1)==iN0 & BCU(:,2)==iDOF),3);
        end
    end
    IndexUPknow(iN1,end)=IndexUPknow(iN1,end)+2; %所有位移连续条件给从节点
end
%施加剪应力正应力相等条件iDOF代表前后
for iCE=1:size(CommonElement,1)
    for iN=1:2
        iNM=Element(CommonElement(iCE,1),2+iN);
        iNS=Element(CommonElement(iCE,3),5-iN);
        ColPosM=2*size(Node,1)-sum(sum(IndexUPknow(:,2:3)))+...
            4*(iNM-1)-sum(sum(IndexUPknow(1:iNM-1,4:7)))+...
            (2-sum(IndexUPknow(iNM,4:5)))*(2-iN);
        ColPosS=2*size(Node,1)-sum(sum(IndexUPknow(:,2:3)))+...
            4*(iNS-1)-sum(sum(IndexUPknow(1:iNS-1,4:7)))+...
            (2-sum(IndexUPknow(iNS,4:5)))*(iN-1);  
        HGL(end+(1:2),[ColPosM+(1:2) ColPosS+(1:2)])=[eye(2) eye(2)];
        F(end+(1:2))=0;
        
        if CommonElement(iCE,1)==Element(find(Element(:,3)==iNM),2)
            E2M=Element(find(Element(:,4)==iNM),2);
        else
            E2M=Element(find(Element(:,3)==iNM),2);
        end
        E2S=CommonElement(find(CommonElement(:,1)==E2M),3);
        bMidNode=0;
        if ~isempty(E2S) & Element(E2S,1)==Element(CommonElement(iCE,3),1)
            bMidNode=1;%共边的中间节点
        end       
        
        if bMidNode==1 
%             |...
%            (IndexUPknow(iNM,end) < 4 & IndexUPknow(iNS,end) < 4)
            IndexUPknow(iNM,end)=IndexUPknow(iNM,end)+1;
            IndexUPknow(iNS,end)=IndexUPknow(iNS,end)+1;
        else
            for iDOF=1:2        %所有的连续条件首先满足方程少的点
%                 以前的程序
%                 if IndexUPknow(iNM,end)<= IndexUPknow(iNS,end)
%                     IndexUPknow(iNM,end)=IndexUPknow(iNM,end)+1; 
%                 else
%                     IndexUPknow(iNS,end)=IndexUPknow(iNS,end)+1;
%                 end
                Rows=[iNM iNS];
                [CMin IMin] = min(IndexUPknow(Rows,end));
                RMin=Rows(IMin);
                IndexUPknow(RMin,end)=IndexUPknow(RMin,end)+1;
            end            
        end
    end 
end
%分析、判断----补充方程
NodeNo=find(IndexUPknow(:,end)==3 | IndexUPknow(:,end)==2);
for iN=1:numel(NodeNo)
    Nodei=NodeNo(iN);
    Eb=find(Element(:,4)==Nodei);
    Ea=find(Element(:,3)==Nodei);
    N1=NE(Eb,:);
    N2=NE(Ea,:);
    N21=[N2 -N1];
    PositonP=2*size(Node,1)-sum(sum(IndexUPknow(:,2:3)))+...
        4*(Nodei-1)-sum(sum(IndexUPknow(1:Nodei-1,4:7)));
    indexP0 = find(IndexUPknow(Nodei,4:7)==0);
    indexP1 = setdiff([1:4],indexP0);
    HGL(end+1,PositonP+(1:numel(indexP0)))=N21(indexP0);
    F(end+1)=-N21(indexP1)*UPknow(Nodei,3+indexP1)';

    if IndexUPknow(Nodei,end)==3
        continue
    end
    
    
    IndexRow=size(HGL,1)+1;
    if  norm(N1-N2)<1e-3  %两单元共线，特殊的角点约束
        MN=(N1+N2)/2;
        N21=[MN(2) -MN(1) -MN(2) MN(1)];
        HGL(IndexRow,PositonP+(1:numel(indexP0)))=N21(indexP0);
        F(IndexRow)=-N21(indexP1)*UPknow(Nodei,3+indexP1)';
    elseif MChau==1%加另一个角点约束            
        G2=2*Mat(Area(Element(Eb,1),2),4)/Gm;
        CChauN1=-G2/LE(Eb)*[-N1(2) N1(1)];
        CChauN2=-G2/LE(Ea)*[-N2(2) N2(1)];
        CChauN0=-CChauN1-CChauN2;
        CChauU=[CChauN1;CChauN0;CChauN2];
        CChauP=[N1 -N2];
        NChau=[Element(Eb,3) Nodei Element(Ea,4)];
        F(IndexRow)=0;
        for iNChau=1:3
            IndexColumn=2*(NChau(iNChau)-1)-sum(sum(IndexUPknow(1:(NChau(iNChau)-1),2:3)));
            for iDof=1:2
                if IndexUPknow(NChau(iNChau),iDof+1)==0 %自由度未知
                    IndexColumn=IndexColumn+1;
                    HGL(IndexRow,IndexColumn)=CChauU(iNChau,iDof);
                else    %自由度已知
                    UiDof=BCU(find(BCU(:,1)==NChau(iNChau) & BCU(:,2)==iDof),3);
                    F(IndexRow)=F(IndexRow)-CChauU(iNChau,iDof)*UiDof;
                end
            end
        end
        HGL(IndexRow,PositonP+(1:4-sum(IndexUPknow(Nodei,4:7))))=CChauP(indexP0);
        F(IndexRow)=F(IndexRow)-CChauP(indexP1)*UPknow(Nodei,3+indexP1)';
    else
        %-------------------虚节点法-------------------begin
        
        HGL(end,:)=[];
        F(end)=[];
        
        ENo=find(Element(:,3)==Nodei);
        Nodej=Element(ENo,4);
        Nodek=size(Node,1)+1;
        Node(end+1,:)=[size(Node,1)+1 sum(Node([Nodei Nodej] ,2:3))/2];
        
        Element(end+1,:)=Element(ENo,:);
        Element(ENo,4)=Nodek;
        Element(end,[2 3])=[size(Element,1) Nodek];
        LE(ENo)=LE(ENo)/2;
        LE(end+1)=LE(ENo);
        NE(end+1,:)=NE(ENo,:);        
        
        np=Nodek;
        iA=Element(ENo,1);
        ElementAi=Element(find(Element(:,1)==iA),:);
        NodeAi=ElementAi(:,3:4);
        disp(sprintf('Area=%d        iN=%d',iA,np))
        
        %求C矩阵 
        n1=ElementAi(find(ElementAi(:,4)==np),3);
        n2=ElementAi(find(ElementAi(:,3)==np),4);
        a1=DirectionAngle(Node(np,2:end),Node(n1,2:end));
        a2=DirectionAngle(Node(np,2:end),Node(n2,2:end));
        if a2<a1
            a2=a2+2*pi;
        end
        da=a2-a1;
        aa=(a1+a2);
        nu2=2*(1-NU);
        Ci=eye(2)-1/2/pi/nu2*[nu2*da+cos(aa)*sin(da) sin(aa)*sin(da);...
                              sin(aa)*sin(da) nu2*da-cos(aa)*sin(da)];
        RIndex=2*np+[-1 0];
        C(RIndex,RIndex)=Ci;
        
        xp=Node(np,2:end)';
        [Hnp Gnp VFnp]=HGFA(np,ElementAi,IntMethod);
        
        Hnp(:,end-1:end)=Hnp(:,end-1:end)+Ci;
        HalfOneN=0.5*[diag(ones(2,1)) diag(ones(2,1))];
        Hnp(:,[Nodei*2+(-1:0) Nodej*2+(-1:0)])=...
            Hnp(:,[Nodei*2+(-1:0) Nodej*2+(-1:0) Nodek*2+(-1:0)])*...
            [diag(ones(4,1));HalfOneN];
        Hnp(:,Nodek*2+(-1:0))=[];

        Gnp=Gnp*Gm;
        Gnp(:,[Nodei*4+(-1:0) Nodej*4+(-3:-2)])=...
            Gnp(:,[Nodei*4+(-1:0) Nodej*4+(-3:-2) Nodek*4+(-3:0)])*...
            [diag(ones(4,1));HalfOneN;HalfOneN];
        Gnp(:,Nodek*4+(-3:0))=[];
        
        %计算
        HGL(end+1:end+2,:)=cat(2,Hnp(:,IndexU0),-Gnp(:,IndexP0));
        HGRnp=cat(2,-Hnp(:,IndexU1),Gnp(:,IndexP1));
        HGR(end+1:end+2,:)=HGRnp;        
        F(end+1:end+2)=HGRnp*UP1+VFnp;
        

        
        Node(end,:)=[];
        Element(end,:)=[];
        Element(ENo,4)=Nodej;
        LE(ENo)=LE(ENo)*2;
        LE(end)=[];
        NE(end,:)=[];
        %-------------------虚节点法-------------------end
        
        %-------------------临单元剪应力约束法-------------------begin
%         Nodej=Element(find(Element(:,3)==Nodei),4);
%         Nodek=Element(find(Element(:,3)==Nodej),4);
%         PositonP=2*size(Node,1)-sum(sum(IndexUPknow(:,2:3)))+...
%             4*(Nodei-1)-sum(sum(IndexUPknow(1:Nodei-1,4:7)));
%         Posi=PositonP;
%         PositonP=2*size(Node,1)-sum(sum(IndexUPknow(:,2:3)))+...
%             4*(Nodej-1)-sum(sum(IndexUPknow(1:Nodej-1,4:7)));
%         Posj=PositonP;
%         PositonP=2*size(Node,1)-sum(sum(IndexUPknow(:,2:3)))+...
%             4*(Nodek-1)-sum(sum(IndexUPknow(1:Nodek-1,4:7)));
%         Posk=PositonP;
%         HGL(end+1,[Posi Posj Posk])=[1 -2 1];
%          F(end+1)=0;
%         HGL(end+1,[Posi Posj Posk]+1)=[1 -2 1];
%          F(end+1)=0;
         %-------------------临单元剪应力约束法-------------------end
    end 
end

%C------------------------------
UP0=HGL\F;
U0=UP0(1:(2*size(Node,1)-sum(sum(IndexUPknow(:,2:3)))));
P0=UP0((2*size(Node,1)-sum(sum(IndexUPknow(:,2:3))))+1:end);
P0=P0*Gm;
Utemp=zeros(size(Node,1),2)';
Utemp(find(IndexUPknow(:,2:3)'==0))=U0;
UPknow(:,2:3)=UPknow(:,2:3)+Utemp';
Ptemp=zeros(size(Node,1),4)';
Ptemp(find(IndexUPknow(:,4:7)'==0))=P0;
UPknow(:,4:7)=UPknow(:,4:7)+Ptemp';
disp('boundary point:');
disp(sprintf('%s\t\t\t\t%s\t\t\t\t\t\t\t%s','Node No','UX','UY'));
disp(sprintf('%d\t%24.12f\t%24.12f\n',UPknow(:,1:3)'));
disp(sprintf('%s\t\t\t%s\t\t\t\t\t%s\t\t\t\t\t%s\t\t\t\t\t%s','Node No','SXb','SYb','SXa','SYa'))
disp(sprintf('%d\t%16.4f\t%16.4f\t%16.4f\t%16.4f\n',UPknow(:,[1 4:7])'))


% 计算内部点的位移
InCount=size(InNode,1);
iNterU=zeros(0,5);
iNterP=zeros(0,6);
for iKP=1:InCount
    xp=(InNode(iKP,1:2))';
    iKPState=-1;           %IN(any area)=1  ON=0  OHTERS=-1
    for iA=1:max(Element(:,1))
        MatiA=Area(iA,2);
        Den=Mat(MatiA,1);
        NU=Mat(MatiA,3);
        SG=Mat(MatiA,4);   
        EA=Element(find(Element(:,1)==iA),3:4);
        LoopNo=1;
        ON=0;
        while size(EA,1)~=0
            
            LoopNode=zeros(0,1);
            [Cmin Imin]=min(Node(EA(:,1),2));
            LoopNode(1:2,1)=EA(Imin,:)';
            EA(Imin,:)=[];
            while LoopNode(1)~=LoopNode(end)
                IRow=find(EA(:,1)==LoopNode(end));
                LoopNode(end+1)=EA(IRow,2);
                EA(IRow,:)=[];
            end
            LoopNode(end)=[];
            [INLoop ONLoop] = inpolygon(xp(1),xp(2),Node(LoopNode,2),Node(LoopNode,3));
            if ONLoop==1
                ON=1;
                break
            end            
            if LoopNo==1 
                if INLoop==0
                    IN=0;
                    break
                else
                    IN=1;
                end
            else
                if INLoop==1
                    IN==0;
                    break
                end
            end
            LoopNo=LoopNo+1;
        end

        ElementAi=Element(find(Element(:,1)==iA),:);
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
% 单元法向与y轴的夹角
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
                if y0<0
                end
                Constp=LE(ElNo)/8/pi/SG/(1-NU)/(x2-x1)^2;
                Constq=-LE(ElNo)/4/pi/(1-NU)/(x2-x1)^2;
                Constr=LE(ElNo)/4/pi/(1-NU)/(x2-x1)^2;
                Consts=LE(ElNo)*SG/2/pi/(1-NU)/(x2-x1)^2;
%                     (7-8*NU)/2
                temp1=-(7-8*NU)/2*(x2*(x2-x1)-(x2^2-x1^2)/2);
                temp2=-(7-8*NU)/2*(-x1*(x2-x1)+(x2^2-x1^2)/2);
                n1=NE(ElNo,1);
                n2=NE(ElNo,2);
                BLOG=0;
% Ui            
                    if abs(y0)<1e-20
                        if abs(x2)<1e-30
                            Treat2=x1^2/2;
                            Treat6=x1^2/2;
                            y2Const2=0;
                            y2Const7=0;
                            yConst6=0;
                            yConst11=0;
                            B1=(4*NU-3)/2*(-1/2*(-x1^2*log(x1^2)+x1^2));
                            B5=(4*NU-3)/2*(-x1*(-x1*log(x1^2)+2*x1)+1/2*(-x1^2*log(x1^2)+x1^2));
                            yConst2=0;
                            yConst7=-x1*pi/2*(-sign(x1));
                            yConst3=0;
                            yConst8=-x1/2*pi/2*(-sign(x1));
                            y3Const4=0; 
                            y3Const9=-x1/2*pi/2*(-sign(x1)); 
                            y2Const5=0;
                            y2Const10=0;
                            limConst6=x1; 
                            limConst11=x1/2*(log(x1^2))-x1; 
                            
                            B2=a^2*Treat2+b^2*y2Const2+2*a*b*yConst6;
                            B3=c^2*Treat2+a^2*y2Const2+2*a*c*yConst6;
                            B4=a*c*Treat2+(a^2+b*c)*yConst6+a*b*y2Const2;  
                            B6=a^2*Treat6+b^2*y2Const7+2*a*b*yConst11;
                            B7=c^2*Treat6+a^2*y2Const7+2*a*c*yConst11;
                            B8=a*c*Treat6+(a^2+b*c)*yConst11+a*b*y2Const7;
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
                        
    %                     Ui
                                                   
                            A1=(1-2*NU)*yConst2;
                            A2=2*(a^2*yConst3+b^2*y3Const4+2*a*b*y2Const5);
                            A3=2*(c^2*yConst3+a^2*y3Const4+2*a*c*y2Const5);
                            A4=2*(a*c*yConst3+a*b*y3Const4+(a^2+b*c)*y2Const5);
                            A5=-(1-2*NU)*((a*n2-c*n1)*limConst6+(b*n2-a*n1)*yConst2);
                            A6=(1-2*NU)*yConst7;
                            A7=2*(a^2*yConst8+b^2*y3Const9+2*a*b*y2Const10);
                            A8=2*(c^2*yConst8+a^2*y3Const9+2*a*c*y2Const10);
                            A9=2*(a*c*yConst8+a*b*y3Const9+(a^2+b*c)*y2Const10);
                            A10=-(1-2*NU)*((a*n2-c*n1)*limConst11+(b*n2-a*n1)*yConst7);
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
                        
                            rou1=a*limConst6+b*yConst2;
                            je1=limConst6-y2Const5;
                            rou2=a^3*je1+3*a^2*b*yConst3+3*a*b^2*y2Const5+b^3*y3Const4;
                            rou3=a*limConst11+b*yConst7;
                            je2=limConst11-y2Const10;
                            rou4=a^3*je2+3*a^2*b*yConst8+3*a*b^2*y2Const10+b^3*y3Const9;
                            rou5=-c*limConst6-a*yConst2;
                            rou6=a^2*c*je1+(a^3+2*a*b*c)*yConst3+(b^2*c+2*a^2*b)*y2Const5+a*b^2*y3Const4;
                            rou7=-c*limConst11-a*yConst7;
                            rou8=a^2*c*je2+(a^3+2*a*b*c)*yConst8+(b^2*c+2*a^2*b)*y2Const10+a*b^2*y3Const9;
                            rou9=a*c^2*je1+(2*a^2*c+b*c^2)*yConst3+(a^3+2*a*b*c)*y2Const5+a^2*b*y3Const4;
                            rou10=a*c^2*je2+(2*a^2*c+b*c^2)*yConst8+(a^3+2*a*b*c)*y2Const10+a^2*b*y3Const9;
                            rou11=c^3*je1+3*c^2*a*yConst3+3*c*a^2*y2Const5+a^3*y3Const4;
                            rou12=c^3*je2+3*c^2*a*yConst8+3*c*a^2*y2Const10+a^3*y3Const9;
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
                            AD=Constr*[D111q1 D211q1 D111q2 D211q2;D112q1 D212q1 D112q2 D212q2;D122q1 D222q1 D122q2 D222q2];

                           
                            ycst1=0;
                            x2ycst1=0;
                            x1ycst1=0;
                            y2cst1=0;
                            x2y2cst2=1/8;
                            x1y2cst2=-1/8;
                            y3cst2=1/8*pi/2*(-sign(x1));
                            y3cst3=0;
                            x2y3cst3=0;
                            x1y3cst3=0;
                            y4cst3=0; 
                            ycst4=3/8*pi/2*(-sign(x1));
                            x2y4h2=3/8;
                            x1y4h2=-3/8;
                            yConst5=-pi/4*(-sign(x1));     
                            yConst10=pi/4*(-sign(x1));
                            y2Const4=1/2;
                            y2Const9=1/2;
                            limConst2=log(x1^2)/2;
                            limConst7=-log(x1^2)/2;
                            limConst3=-1/2+log(x1^2)/2;
                            limConst8=-1/2-log(x1^2)/2;
                           
                           ywer1=a*yConst5+b*y2Const4;
                           ywer2=a*yConst10+b*y2Const9;
                           ywer3=c*yConst5+a*y2Const4;
                           ywer4=c*yConst10+a*y2Const9;
                           ymj1=a^3*(x2ycst1-ycst4)+3*a^2*b*(x2y2cst2-y2cst1)+3*a*b^2*(x2y3cst3-y3cst2)+b^3*(x2y4h2-y4cst3);
                           ymj2=a^3*(ycst4-x1ycst1)+3*a^2*b*(y2cst1-x1y2cst2)+3*a*b^2*(y3cst2-x1y3cst3)+b^3*(y4cst3-x1y4h2);
                           ymj3=a^2*c*(x2ycst1-ycst4)+(a^3+2*a*b*c)*(x2y2cst2-y2cst1)+(b^2*c+2*a^2*b)*(x2y3cst3-y3cst2)+b^2*a*(x2y4h2-y4cst3);
                           ymj4=a^2*c*(ycst4-x1ycst1)+(a^3+2*a*b*c)*(y2cst1-x1y2cst2)+(b^2*c+2*a^2*b)*(y3cst2-x1y3cst3)+b^2*a*(y4cst3-x1y4h2);
                           ymj5=a*c^2*(x2ycst1-ycst4)+(2*a^2*c+b*c^2)*(x2y2cst2-y2cst1)+(a^3+2*a*b*c)*(x2y3cst3-y3cst2)+a^2*b*(x2y4h2-y4cst3);
                           ymj6=a*c^2*(ycst4-x1ycst1)+(2*a^2*c+b*c^2)*(y2cst1-x1y2cst2)+(a^3+2*a*b*c)*(y3cst2-x1y3cst3)+a^2*b*(y4cst3-x1y4h2);
                           ymj7=c^3*(x2ycst1-ycst4)+3*c^2*a*(x2y2cst2-y2cst1)+3*c*a^2*(x2y3cst3-y3cst2)+a^3*(x2y4h2-y4cst3);
                           ymj8=c^3*(ycst4-x1ycst1)+3*c^2*a*(y2cst1-x1y2cst2)+3*c*a^2*(y3cst2-x1y3cst3)+a^3*(y4cst3-x1y4h2);
                           yme1=ywer1-4*ymj1;
                           me2=a^2*limConst3+2*a*b*yConst5+b^2*y2Const4;
                           yme3=ywer2-4*ymj2;
                           me4=a^2*limConst8+2*a*b*yConst10+b^2*y2Const9;
                           yme5=(1-2*NU)*ywer3-4*ymj3;
                           me6=a*c*limConst3+(a^2+b*c)*yConst5+a*b*y2Const4;
                           yme7=(1-2*NU)*ywer4-4*ymj4;
                           me8=a*c*limConst8+(a^2+b*c)*yConst10+a*b*y2Const9;
                           yme9=NU*ywer3-4*ymj3;
                           yme10=NU*ywer4-4*ymj4;
                           yme11=NU*ywer1-4*ymj5;
                           me12=c^2*limConst3+2*a*c*yConst5+a^2*y2Const4;
                           yme13=NU*ywer2-4*ymj6;
                           me14=c^2*limConst8+2*a*c*yConst10+a^2*y2Const9;
                           yme15=(1-2*NU)*ywer1-4*ymj5;
                           yme16=(1-2*NU)*ywer2-4*ymj6;
                           yme17=ywer3-4*ymj7;
                           yme18=ywer4-4*ymj8;
                           S111q1=2*yme1+2*n1*me2+n1*limConst2;
                           S111q2=2*yme3+2*n1*me4+n1*limConst7;
                           S211q1=2*yme5+4*NU*n1*me6+(1-2*NU)*2*n2*me2-(1-4*NU)*n2*limConst2;
                           S211q2=2*yme7+4*NU*n1*me8+(1-2*NU)*2*n2*me4-(1-4*NU)*n2*limConst7;
                           S112q1=2*yme9+(2-2*NU)*n1*me6+2*NU*n2*me2+(1-2*NU)*n2*limConst2;
                           S112q2=2*yme10+(2-2*NU)*n1*me8+2*NU*n2*me4+(1-2*NU)*n2*limConst7;
                           S212q1=2*yme11+2*n1*NU*me12+(2-2*NU)*n2*me6+(1-2*NU)*n1*limConst2;
                           S212q2=2*yme13+2*n1*NU*me14+(2-2*NU)*n2*me8+(1-2*NU)*n1*limConst7;
                           S122q1=2*yme15+4*NU*n2*me6+(1-2*NU)*2*n1*me12-(1-4*NU)*n1*limConst2;
                           S122q2=2*yme16+4*NU*n2*me8+(1-2*NU)*2*n1*me14-(1-4*NU)*n1*limConst7;
                           S222q1=2*yme17+2*n2*me12+n2*limConst2;
                           S222q2=2*yme18+2*n2*me14+n2*limConst7;
                           AS=Consts*[S111q1 S211q1 S111q2 S211q2;S112q1 S212q1 S112q2 S212q2;S122q1 S222q1 S122q2 S222q2];
                        elseif abs(x1)<1e-30
                           Treat2=x2^2/2;
                            Treat6=x2^2/2;
                            y2Const2=0;
                            y2Const7=0;
                            yConst6=0;
                            yConst11=0;
                            B1=(4*NU-3)/2*(x2*(x2*log(x2^2)-2*x2)-1/2*(x2^2*log(x2^2)-x2^2));
                            B5=(4*NU-3)/2*(1/2*(x2^2*log(x2^2)-x2^2));
                            % yConst2=x2*(atan(x2/y0)-atan(x1/y0));
                            yConst2=x2*pi/2*sign(x2);
                            yConst7=0;
                            yConst3=x2/2*pi/2*sign(x2);
                            yConst8=0;
                            y3Const4=x2/2*pi/2*sign(x2); 
                            y3Const9=0; 
                            y2Const5=0;
                            y2Const10=0;
                            limConst6=x2/2*log(x2^2)-x2; 
                            limConst11=x2; 
                           
                            B2=a^2*Treat2+b^2*y2Const2+2*a*b*yConst6;
                            B3=c^2*Treat2+a^2*y2Const2+2*a*c*yConst6;
                            B4=a*c*Treat2+(a^2+b*c)*yConst6+a*b*y2Const2;                           
                            B6=a^2*Treat6+b^2*y2Const7+2*a*b*yConst11;
                            B7=c^2*Treat6+a^2*y2Const7+2*a*c*yConst11;
                            B8=a*c*Treat6+(a^2+b*c)*yConst11+a*b*y2Const7;
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
                        
    %                     Ui                       
                            A1=(1-2*NU)*yConst2;
                            A2=2*(a^2*yConst3+b^2*y3Const4+2*a*b*y2Const5);
                            A3=2*(c^2*yConst3+a^2*y3Const4+2*a*c*y2Const5);
                            A4=2*(a*c*yConst3+a*b*y3Const4+(a^2+b*c)*y2Const5);
                            A5=-(1-2*NU)*((a*n2-c*n1)*limConst6+(b*n2-a*n1)*yConst2);
                            A6=(1-2*NU)*yConst7;
                            A7=2*(a^2*yConst8+b^2*y3Const9+2*a*b*y2Const10);
                            A8=2*(c^2*yConst8+a^2*y3Const9+2*a*c*y2Const10);
                            A9=2*(a*c*yConst8+a*b*y3Const9+(a^2+b*c)*y2Const10);
                            A10=-(1-2*NU)*((a*n2-c*n1)*limConst11+(b*n2-a*n1)*yConst7);
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
                        
                            rou1=a*limConst6+b*yConst2;
                            je1=limConst6-y2Const5;
                            rou2=a^3*je1+3*a^2*b*yConst3+3*a*b^2*y2Const5+b^3*y3Const4;
                            rou3=a*limConst11+b*yConst7;
                            je2=limConst11-y2Const10;
                            rou4=a^3*je2+3*a^2*b*yConst8+3*a*b^2*y2Const10+b^3*y3Const9;
                            rou5=-c*limConst6-a*yConst2;
                            rou6=a^2*c*je1+(a^3+2*a*b*c)*yConst3+(b^2*c+2*a^2*b)*y2Const5+a*b^2*y3Const4;
                            rou7=-c*limConst11-a*yConst7;
                            rou8=a^2*c*je2+(a^3+2*a*b*c)*yConst8+(b^2*c+2*a^2*b)*y2Const10+a*b^2*y3Const9;
                            rou9=a*c^2*je1+(2*a^2*c+b*c^2)*yConst3+(a^3+2*a*b*c)*y2Const5+a^2*b*y3Const4;
                            rou10=a*c^2*je2+(2*a^2*c+b*c^2)*yConst8+(a^3+2*a*b*c)*y2Const10+a^2*b*y3Const9;
                            rou11=c^3*je1+3*c^2*a*yConst3+3*c*a^2*y2Const5+a^3*y3Const4;
                            rou12=c^3*je2+3*c^2*a*yConst8+3*c*a^2*y2Const10+a^3*y3Const9;
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
                            AD=Constr*[D111q1 D211q1 D111q2 D211q2;D112q1 D212q1 D112q2 D212q2;D122q1 D222q1 D122q2 D222q2];

              
                            ycst1=0;
                            x2ycst1=0;
                            x1ycst1=0;
                            y2cst1=0;
                            x2y2cst2=1/8;
                            x1y2cst2=-1/8;
                            y3cst2=1/8*pi/2*sign(x2);
                            y3cst3=0;
                            x2y3cst3=0;
                            x1y3cst3=0;
                            y4cst3=0; 
                            ycst4=3/8*pi/2*sign(x2);
                            x2y4h2=3/8;
                            x1y4h2=-3/8;
                            yConst5=-pi/4*sign(x2);     
                            yConst10=pi/4*sign(x2);
                            y2Const4=1/2;
                            y2Const9=1/2;
                            limConst2=-log(x2^2)/2;
                            limConst7=log(x2^2)/2;
                            limConst3=-1/2-log(x2^2)/2;
                            limConst8=-1/2+log(x2^2)/2;
                            
                           ywer1=a*yConst5+b*y2Const4;
                           ywer2=a*yConst10+b*y2Const9;
                           ywer3=c*yConst5+a*y2Const4;
                           ywer4=c*yConst10+a*y2Const9;
                           ymj1=a^3*(x2ycst1-ycst4)+3*a^2*b*(x2y2cst2-y2cst1)+3*a*b^2*(x2y3cst3-y3cst2)+b^3*(x2y4h2-y4cst3);
                           ymj2=a^3*(ycst4-x1ycst1)+3*a^2*b*(y2cst1-x1y2cst2)+3*a*b^2*(y3cst2-x1y3cst3)+b^3*(y4cst3-x1y4h2);
                           ymj3=a^2*c*(x2ycst1-ycst4)+(a^3+2*a*b*c)*(x2y2cst2-y2cst1)+(b^2*c+2*a^2*b)*(x2y3cst3-y3cst2)+b^2*a*(x2y4h2-y4cst3);
                           ymj4=a^2*c*(ycst4-x1ycst1)+(a^3+2*a*b*c)*(y2cst1-x1y2cst2)+(b^2*c+2*a^2*b)*(y3cst2-x1y3cst3)+b^2*a*(y4cst3-x1y4h2);
                           ymj5=a*c^2*(x2ycst1-ycst4)+(2*a^2*c+b*c^2)*(x2y2cst2-y2cst1)+(a^3+2*a*b*c)*(x2y3cst3-y3cst2)+a^2*b*(x2y4h2-y4cst3);
                           ymj6=a*c^2*(ycst4-x1ycst1)+(2*a^2*c+b*c^2)*(y2cst1-x1y2cst2)+(a^3+2*a*b*c)*(y3cst2-x1y3cst3)+a^2*b*(y4cst3-x1y4h2);
                           ymj7=c^3*(x2ycst1-ycst4)+3*c^2*a*(x2y2cst2-y2cst1)+3*c*a^2*(x2y3cst3-y3cst2)+a^3*(x2y4h2-y4cst3);
                           ymj8=c^3*(ycst4-x1ycst1)+3*c^2*a*(y2cst1-x1y2cst2)+3*c*a^2*(y3cst2-x1y3cst3)+a^3*(y4cst3-x1y4h2);
                           yme1=ywer1-4*ymj1;
                           me2=a^2*limConst3+2*a*b*yConst5+b^2*y2Const4;
                           yme3=ywer2-4*ymj2;
                           me4=a^2*limConst8+2*a*b*yConst10+b^2*y2Const9;
                           yme5=(1-2*NU)*ywer3-4*ymj3;
                           me6=a*c*limConst3+(a^2+b*c)*yConst5+a*b*y2Const4;
                           yme7=(1-2*NU)*ywer4-4*ymj4;
                           me8=a*c*limConst8+(a^2+b*c)*yConst10+a*b*y2Const9;
                           yme9=NU*ywer3-4*ymj3;
                           yme10=NU*ywer4-4*ymj4;
                           yme11=NU*ywer1-4*ymj5;
                           me12=c^2*limConst3+2*a*c*yConst5+a^2*y2Const4;
                           yme13=NU*ywer2-4*ymj6;
                           me14=c^2*limConst8+2*a*c*yConst10+a^2*y2Const9;
                           yme15=(1-2*NU)*ywer1-4*ymj5;
                           yme16=(1-2*NU)*ywer2-4*ymj6;
                           yme17=ywer3-4*ymj7;
                           yme18=ywer4-4*ymj8;
                           S111q1=2*yme1+2*n1*me2+n1*limConst2;
                           S111q2=2*yme3+2*n1*me4+n1*limConst7;
%                            S211q1=+(1-2*NU)*2*n2*me2-(1-4*NU)*n2*limConst2;
                           S211q1=2*yme5+4*NU*n1*me6+(1-2*NU)*2*n2*me2-(1-4*NU)*n2*limConst2;
                           S211q2=2*yme7+4*NU*n1*me8+(1-2*NU)*2*n2*me4-(1-4*NU)*n2*limConst7;
                           S112q1=2*yme9+(2-2*NU)*n1*me6+2*NU*n2*me2+(1-2*NU)*n2*limConst2;
                           S112q2=2*yme10+(2-2*NU)*n1*me8+2*NU*n2*me4+(1-2*NU)*n2*limConst7;
                           S212q1=2*yme11+2*n1*NU*me12+(2-2*NU)*n2*me6+(1-2*NU)*n1*limConst2;
                           S212q2=2*yme13+2*n1*NU*me14+(2-2*NU)*n2*me8+(1-2*NU)*n1*limConst7;
                           S122q1=2*yme15+4*NU*n2*me6+(1-2*NU)*2*n1*me12-(1-4*NU)*n1*limConst2;
                           S122q2=2*yme16+4*NU*n2*me8+(1-2*NU)*2*n1*me14-(1-4*NU)*n1*limConst7;
                           S222q1=2*yme17+2*n2*me12+n2*limConst2;
                           S222q2=2*yme18+2*n2*me14+n2*limConst7;
                           AS=Consts*[S111q1 S211q1 S111q2 S211q2;S112q1 S212q1 S112q2 S212q2;S122q1 S222q1 S122q2 S222q2];
                            
%                         ----------------------y0逼近于0的情况-------------------
                        else
                            Treat2=x2*(x2-x1)-(x2^2-x1^2)/2;
                            Treat6=-x1*(x2-x1)+(x2^2-x1^2)/2;
                            y2Const2=0;
                            y2Const7=0;
                            yConst6=0;
                            yConst11=0;
                            B1=(4*NU-3)/2*(x2*(x2*log(x2^2)-x1*log(x1^2)-2*(x2-x1))-...
                                1/2*(x2^2*log(x2^2)-x1^2*log(x1^2)-(x2^2-x1^2)));
                            B5=(4*NU-3)/2*(-x1*(x2*log(x2^2)-x1*log(x1^2)-2*(x2-x1))+...
                               1/2*(x2^2*log(x2^2)-x1^2*log(x1^2)-(x2^2-x1^2)));
                            % yConst2=x2*(atan(x2/y0)-atan(x1/y0));
                            yConst2=x2*pi/2*(sign(x2)-sign(x1));
                            yConst7=-x1*pi/2*(sign(x2)-sign(x1));
                            yConst3=x2/2*pi/2*(sign(x2)-sign(x1));
                            yConst8=-x1/2*pi/2*(sign(x2)-sign(x1));
                            y3Const4=x2/2*pi/2*(sign(x2)-sign(x1)); 
                            y3Const9=-x1/2*pi/2*(sign(x2)-sign(x1)); 
                            y2Const5=0;
                            y2Const10=0;
                            limConst6=x2/2*log((x2^2)/(x1^2))-(x2-x1); 
                            limConst11=-x1/2*log((x2^2)/(x1^2))+(x2-x1); 
                            
                            B2=a^2*Treat2+b^2*y2Const2+2*a*b*yConst6;
                            B3=c^2*Treat2+a^2*y2Const2+2*a*c*yConst6;
                            B4=a*c*Treat2+(a^2+b*c)*yConst6+a*b*y2Const2;  
                            B6=a^2*Treat6+b^2*y2Const7+2*a*b*yConst11;
                            B7=c^2*Treat6+a^2*y2Const7+2*a*c*yConst11;
                            B8=a*c*Treat6+(a^2+b*c)*yConst11+a*b*y2Const7;
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
                        
    %                     Ui
                            
                        
                            A1=(1-2*NU)*yConst2;
                            A2=2*(a^2*yConst3+b^2*y3Const4+2*a*b*y2Const5);
                            A3=2*(c^2*yConst3+a^2*y3Const4+2*a*c*y2Const5);
                            A4=2*(a*c*yConst3+a*b*y3Const4+(a^2+b*c)*y2Const5);
                            A5=-(1-2*NU)*((a*n2-c*n1)*limConst6+(b*n2-a*n1)*yConst2);
                            A6=(1-2*NU)*yConst7;
                            A7=2*(a^2*yConst8+b^2*y3Const9+2*a*b*y2Const10);
                            A8=2*(c^2*yConst8+a^2*y3Const9+2*a*c*y2Const10);
                            A9=2*(a*c*yConst8+a*b*y3Const9+(a^2+b*c)*y2Const10);
                            A10=-(1-2*NU)*((a*n2-c*n1)*limConst11+(b*n2-a*n1)*yConst7);
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
                        
                            rou1=a*limConst6+b*yConst2;
                            je1=limConst6-y2Const5;
                            rou2=a^3*je1+3*a^2*b*yConst3+3*a*b^2*y2Const5+b^3*y3Const4;
                            rou3=a*limConst11+b*yConst7;
                            je2=limConst11-y2Const10;
                            rou4=a^3*je2+3*a^2*b*yConst8+3*a*b^2*y2Const10+b^3*y3Const9;
                            rou5=-c*limConst6-a*yConst2;
                            rou6=a^2*c*je1+(a^3+2*a*b*c)*yConst3+(b^2*c+2*a^2*b)*y2Const5+a*b^2*y3Const4;
                            rou7=-c*limConst11-a*yConst7;
                            rou8=a^2*c*je2+(a^3+2*a*b*c)*yConst8+(b^2*c+2*a^2*b)*y2Const10+a*b^2*y3Const9;
                            rou9=a*c^2*je1+(2*a^2*c+b*c^2)*yConst3+(a^3+2*a*b*c)*y2Const5+a^2*b*y3Const4;
                            rou10=a*c^2*je2+(2*a^2*c+b*c^2)*yConst8+(a^3+2*a*b*c)*y2Const10+a^2*b*y3Const9;
                            rou11=c^3*je1+3*c^2*a*yConst3+3*c*a^2*y2Const5+a^3*y3Const4;
                            rou12=c^3*je2+3*c^2*a*yConst8+3*c*a^2*y2Const10+a^3*y3Const9;
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
                            AD=Constr*[D111q1 D211q1 D111q2 D211q2;D112q1 D212q1 D112q2 D212q2;D122q1 D222q1 D122q2 D222q2];

                            ycst1=0;
                            x2ycst1=0;
                            x1ycst1=0;
                            y2cst1=0;
                            x2y2cst2=1/8*(1-x2/x1);
                            x1y2cst2=1/8*(x1/x2-1);
                            y3cst2=1/8*pi/2*(sign(x2)-sign(x1));
                            y3cst3=0;
                            x2y3cst3=0;
                            x1y3cst3=0;
                            y4cst3=0; 
                            ycst4=3/8*pi/2*(sign(x2)-sign(x1));
                            x2y4h2=3/8*(1-x2/x1);
                            x1y4h2=3/8*(x1/x2-1);
                            yConst5=-pi/4*(sign(x2)-sign(x1));     
                            yConst10=pi/4*(sign(x2)-sign(x1));
                            y2Const4=x2/2*(1/x2-1/x1);
                            y2Const9=-x1/2*(1/x2-1/x1);
                            limConst2=-(log(x2^2)-log(x1^2))/2;
                            limConst7=(log(x2^2)-log(x1^2))/2;
                            limConst3=-x2/2*(1/x2-1/x1)-log((x2^2)/(x1^2))/2;
                            limConst8=x1/2*(1/x2-1/x1)+log((x2^2)/(x1^2))/2;
                            ywer1=a*yConst5+b*y2Const4;
                           ywer2=a*yConst10+b*y2Const9;
                           ywer3=c*yConst5+a*y2Const4;
                           ywer4=c*yConst10+a*y2Const9;
                           ymj1=a^3*(x2ycst1-ycst4)+3*a^2*b*(x2y2cst2-y2cst1)+3*a*b^2*(x2y3cst3-y3cst2)+b^3*(x2y4h2-y4cst3);
                           ymj2=a^3*(ycst4-x1ycst1)+3*a^2*b*(y2cst1-x1y2cst2)+3*a*b^2*(y3cst2-x1y3cst3)+b^3*(y4cst3-x1y4h2);
                           ymj3=a^2*c*(x2ycst1-ycst4)+(a^3+2*a*b*c)*(x2y2cst2-y2cst1)+(b^2*c+2*a^2*b)*(x2y3cst3-y3cst2)+b^2*a*(x2y4h2-y4cst3);
                           ymj4=a^2*c*(ycst4-x1ycst1)+(a^3+2*a*b*c)*(y2cst1-x1y2cst2)+(b^2*c+2*a^2*b)*(y3cst2-x1y3cst3)+b^2*a*(y4cst3-x1y4h2);
                           ymj5=a*c^2*(x2ycst1-ycst4)+(2*a^2*c+b*c^2)*(x2y2cst2-y2cst1)+(a^3+2*a*b*c)*(x2y3cst3-y3cst2)+a^2*b*(x2y4h2-y4cst3);
                           ymj6=a*c^2*(ycst4-x1ycst1)+(2*a^2*c+b*c^2)*(y2cst1-x1y2cst2)+(a^3+2*a*b*c)*(y3cst2-x1y3cst3)+a^2*b*(y4cst3-x1y4h2);
                           ymj7=c^3*(x2ycst1-ycst4)+3*c^2*a*(x2y2cst2-y2cst1)+3*c*a^2*(x2y3cst3-y3cst2)+a^3*(x2y4h2-y4cst3);
                           ymj8=c^3*(ycst4-x1ycst1)+3*c^2*a*(y2cst1-x1y2cst2)+3*c*a^2*(y3cst2-x1y3cst3)+a^3*(y4cst3-x1y4h2);
                           yme1=ywer1-4*ymj1;
                           me2=a^2*limConst3+2*a*b*yConst5+b^2*y2Const4;
                           yme3=ywer2-4*ymj2;
                           me4=a^2*limConst8+2*a*b*yConst10+b^2*y2Const9;
                           yme5=(1-2*NU)*ywer3-4*ymj3;
                           me6=a*c*limConst3+(a^2+b*c)*yConst5+a*b*y2Const4;
                           yme7=(1-2*NU)*ywer4-4*ymj4;
                           me8=a*c*limConst8+(a^2+b*c)*yConst10+a*b*y2Const9;
                           yme9=NU*ywer3-4*ymj3;
                           yme10=NU*ywer4-4*ymj4;
                           yme11=NU*ywer1-4*ymj5;
                           me12=c^2*limConst3+2*a*c*yConst5+a^2*y2Const4;
                           yme13=NU*ywer2-4*ymj6;
                           me14=c^2*limConst8+2*a*c*yConst10+a^2*y2Const9;
                           yme15=(1-2*NU)*ywer1-4*ymj5;
                           yme16=(1-2*NU)*ywer2-4*ymj6;
                           yme17=ywer3-4*ymj7;
                           yme18=ywer4-4*ymj8;
                           S111q1=2*yme1+2*n1*me2+n1*limConst2;
                           S111q2=2*yme3+2*n1*me4+n1*limConst7;
                           S211q1=2*yme5+4*NU*n1*me6+(1-2*NU)*2*n2*me2-(1-4*NU)*n2*limConst2;
                           S211q2=2*yme7+4*NU*n1*me8+(1-2*NU)*2*n2*me4-(1-4*NU)*n2*limConst7;
                           S112q1=2*yme9+(2-2*NU)*n1*me6+2*NU*n2*me2+(1-2*NU)*n2*limConst2;
                           S112q2=2*yme10+(2-2*NU)*n1*me8+2*NU*n2*me4+(1-2*NU)*n2*limConst7;
                           S212q1=2*yme11+2*n1*NU*me12+(2-2*NU)*n2*me6+(1-2*NU)*n1*limConst2;
                           S212q2=2*yme13+2*n1*NU*me14+(2-2*NU)*n2*me8+(1-2*NU)*n1*limConst7;
                           S122q1=2*yme15+4*NU*n2*me6+(1-2*NU)*2*n1*me12-(1-4*NU)*n1*limConst2;
                           S122q2=2*yme16+4*NU*n2*me8+(1-2*NU)*2*n1*me14-(1-4*NU)*n1*limConst7;
                           S222q1=2*yme17+2*n2*me12+n2*limConst2;
                           S222q2=2*yme18+2*n2*me14+n2*limConst7;
                          AS=Consts*[S111q1 S211q1 S111q2 S211q2;S112q1 S212q1 S112q2 S212q2;S122q1 S222q1 S122q2 S222q2]; 
                        end
%                         ---------------------内部点的情况
                    else
                        Treat2=x2*((x2-x1)-y0*(atan(x2/y0)-atan(x1/y0)))-((x2^2-x1^2)-y0^2*log((x2^2+y0^2)/(x1^2+y0^2)))/2;
                        Treat6=-x1*((x2-x1)-y0*(atan(x2/y0)-atan(x1/y0)))+((x2^2-x1^2)/2-y0^2/2*log((x2^2+y0^2)/(x1^2+y0^2)));
                        y2Const2=x2*y0*(atan(x2/y0)-atan(x1/y0))-y0^2*(log(x2^2+y0^2)-log(x1^2+y0^2))/2;
                        y2Const7=-x1*y0*(atan(x2/y0)-atan(x1/y0))+y0^2*(log(x2^2+y0^2)-log(x1^2+y0^2))/2;
                        yConst6=x2/2*y0*log((x2^2+y0^2)/(x1^2+y0^2))-y0*((x2-x1)-y0*(atan(x2/y0)-atan(x1/y0)));
                        yConst11=-x1/2*y0*log((x2^2+y0^2)/(x1^2+y0^2))...
                                 +y0*((x2-x1)-y0*(atan(x2/y0)-atan(x1/y0))); 
                        B1=(4*NU-3)/2*(x2*(x2*log(x2^2+y0^2)-x1*log(x1^2+y0^2)...
                            -2*((x2-x1)-y0*(atan(x2/y0)-atan(x1/y0))))-...
                            1/2*((x2^2+y0^2)*log(x2^2+y0^2)-(x1^2+y0^2)*log(x1^2+y0^2)-(x2^2-x1^2)));
                        B5=(4*NU-3)/2*(-x1*(x2*log(x2^2+y0^2)-x1*log(x1^2+y0^2)...
                           -2*((x2-x1)-y0*(atan(x2/y0)-atan(x1/y0))))+...
                           1/2*((x2^2+y0^2)*log(x2^2+y0^2)-(x1^2+y0^2)*log(x1^2+y0^2)-(x2^2-x1^2)));
                        yConst2=x2*(atan(x2/y0)-atan(x1/y0))-y0*(log(x2^2+y0^2)-log(x1^2+y0^2))/2;
                        yConst7=-x1*(atan(x2/y0)-atan(x1/y0))+y0*(log(x2^2+y0^2)-log(x1^2+y0^2))/2;
                        yConst3=x2/2*((atan(x2/y0)-atan(x1/y0))-y0*(x2/(x2^2+y0^2)-x1/(x1^2+y0^2)))...
                                -(y0*log((x2^2+y0^2)/(x1^2+y0^2))+y0^3*(1/(x2^2+y0^2)-1/(x1^2+y0^2)))/2;
                        yConst8=-x1/2*((atan(x2/y0)-atan(x1/y0))-y0*(x2/(x2^2+y0^2)-x1/(x1^2+y0^2)))...
                                     +(y0*log((x2^2+y0^2)/(x1^2+y0^2))+y0^3*(1/(x2^2+y0^2)-1/(x1^2+y0^2)))/2;
                        y3Const4=x2/2*(y0*(x2/(x2^2+y0^2)-x1/(x1^2+y0^2))+(atan(x2/y0)-atan(x1/y0)))...
                                      +y0^3*(1/(x2^2+y0^2)-1/(x1^2+y0^2))/2;  
                        y3Const9=-x1/2*(y0*(x2/(x2^2+y0^2)-x1/(x1^2+y0^2))+(atan(x2/y0)-atan(x1/y0)))...
                                        -y0^3*(1/(x2^2+y0^2)-1/(x1^2+y0^2))/2;
                        y2Const5=(y0^2*(x2-x1)/(x1^2+y0^2)-y0*(atan(x2/y0)-atan(x1/y0)))/2;
                        y2Const10=(y0^2*(x1-x2)/(x2^2+y0^2)+y0*(atan(x2/y0)-atan(x1/y0)))/2;
                        limConst6=x2/2*log((x2^2+y0^2)/(x1^2+y0^2))-((x2-x1)-y0*(atan(x2/y0)-atan(x1/y0))); 
                        limConst11=-x1/2*log((x2^2+y0^2)/(x1^2+y0^2))+(x2-x1)-y0*(atan(x2/y0)-atan(x1/y0)); 


                         B2=a^2*Treat2+b^2*y2Const2+2*a*b*yConst6;
                         B3=c^2*Treat2+a^2*y2Const2+2*a*c*yConst6;
                         B4=a*c*Treat2+(a^2+b*c)*yConst6+a*b*y2Const2;   
                         B6=a^2*Treat6+b^2*y2Const7+2*a*b*yConst11;
                         B7=c^2*Treat6+a^2*y2Const7+2*a*c*yConst11;
                         B8=a*c*Treat6+(a^2+b*c)*yConst11+a*b*y2Const7;
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
                        
                        A1=(1-2*NU)*yConst2;
                        A2=2*(a^2*yConst3+b^2*y3Const4+2*a*b*y2Const5);
                        A3=2*(c^2*yConst3+a^2*y3Const4+2*a*c*y2Const5);
                        A4=2*(a*c*yConst3+a*b*y3Const4+(a^2+b*c)*y2Const5);
                        A5=-(1-2*NU)*((a*n2-c*n1)*limConst6+(b*n2-a*n1)*yConst2);
                        A6=(1-2*NU)*yConst7;
                        A7=2*(a^2*yConst8+b^2*y3Const9+2*a*b*y2Const10);
                        A8=2*(c^2*yConst8+a^2*y3Const9+2*a*c*y2Const10);
                        A9=2*(a*c*yConst8+a*b*y3Const9+(a^2+b*c)*y2Const10);
                        A10=-(1-2*NU)*((a*n2-c*n1)*limConst11+(b*n2-a*n1)*yConst7);
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
                    
                        rou1=a*limConst6+b*yConst2;
                        je1=limConst6-y2Const5;
                        rou2=a^3*je1+3*a^2*b*yConst3+3*a*b^2*y2Const5+b^3*y3Const4;
                        rou3=a*limConst11+b*yConst7;
                        je2=limConst11-y2Const10;
                        rou4=a^3*je2+3*a^2*b*yConst8+3*a*b^2*y2Const10+b^3*y3Const9;
                        rou5=-c*limConst6-a*yConst2;
                        rou6=a^2*c*je1+(a^3+2*a*b*c)*yConst3+(b^2*c+2*a^2*b)*y2Const5+a*b^2*y3Const4;
                        rou7=-c*limConst11-a*yConst7;
                        rou8=a^2*c*je2+(a^3+2*a*b*c)*yConst8+(b^2*c+2*a^2*b)*y2Const10+a*b^2*y3Const9;
                        rou9=a*c^2*je1+(2*a^2*c+b*c^2)*yConst3+(a^3+2*a*b*c)*y2Const5+a^2*b*y3Const4;
                        rou10=a*c^2*je2+(2*a^2*c+b*c^2)*yConst8+(a^3+2*a*b*c)*y2Const10+a^2*b*y3Const9;
                        rou11=c^3*je1+3*c^2*a*yConst3+3*c*a^2*y2Const5+a^3*y3Const4;
                        rou12=c^3*je2+3*c^2*a*yConst8+3*c*a^2*y2Const10+a^3*y3Const9;
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
                        AD=Constr*[D111q1 D211q1 D111q2 D211q2;D112q1 D212q1 D112q2 D212q2;D122q1 D222q1 D122q2 D222q2];
%                 ---ASi-
              
                        ycst1=1/2*(y0*(x1^2+y0^2)-y0*(x2^2+y0^2)+y0^3/2*(1/(x2^2+y0^2)^2-1/(x1^2+y0^2)^2));
                        x2ycst1=1/2*(x2*y0/(x1^2+y0^2)-x2*y0/(x2^2+y0^2)+y0^3/2*(x2/(x2^2+y0^2)^2-x2/(x1^2+y0^2)^2));
                        x1ycst1=1/2*(x1*y0/(x1^2+y0^2)-x1*y0/(x2^2+y0^2)+y0^3/2*(x1/(x2^2+y0^2)^2-x1/(x1^2+y0^2)^2));
                        y2cst1=1/2*(y0^2/(x1^2+y0^2)-y0^2/(x2^2+y0^2)+y0^4/2*(1/(x2^2+y0^2)^2-1/(x1^2+y0^2)^2));
                        x2y2cst2=1/8*(x2^2/(x2^2+y0^2)-x1*x2/(x1^2+y0^2)+1/y0*x2*(atan(x2/y0)-atan(x1/y0)))...
                                 -1/4*(x2^2*y0^2/(x2^2+y0^2)^2-x1*x2*y0^2/(x1^2+y0^2)^2);
                        x1y2cst2=1/8*(x2*x1/(x2^2+y0^2)-x1^2/(x1^2+y0^2)+1/y0*x1*(atan(x2/y0)-atan(x1/y0)))...
                                  -1/4*(x2*x1*y0^2/(x2^2+y0^2)^2-x1^2*y0^2/(x1^2+y0^2)^2);
                        y3cst2=1/8*(x2*y0/(x2^2+y0^2)-x1*y0/(x1^2+y0^2)+(atan(x2/y0)-atan(x1/y0)))...
                                 -1/4*(x2*y0^3/(x2^2+y0^2)^2-x1*y0^3/(x1^2+y0^2)^2); 
                        y3cst3=-1/4*y0^3*(1/(x2^2+y0^2)^2-1/(x1^2+y0^2)^2);
                        x2y3cst3=-1/4*x2*y0^3*(1/(x2^2+y0^2)^2-1/(x1^2+y0^2)^2); 
                        x1y3cst3=-1/4*x1*y0^3*(1/(x2^2+y0^2)^2-1/(x1^2+y0^2)^2);
                        y4cst3=-1/4*y0^4*(1/(x2^2+y0^2)^2-1/(x1^2+y0^2)^2); 
                        ycst4=3/8*(atan(x2/y0)-atan(x1/y0))-5/8*(x2*y0/(x2^2+y0^2)-x1*y0/(x1^2+y0^2))...
                               +y0^3/4*(x2/(x2^2+y0^2)^2-x1/(x1^2+y0^2)^2);
                        x2y4h2=1/4*(x2^2*y0^2/(x2^2+y0^2)^2-x1*x2*y0^2/(x1^2+y0^2)^2+...
                               3/2*(x2^2/(x2^2+y0^2)-x1*x2/(x1^2+y0^2)+1/y0*x2*(atan(x2/y0)-atan(x1/y0))));
                        x1y4h2=1/4*(x2*x1*y0^2/(x2^2+y0^2)^2-x1^2*y0^2/(x1^2+y0^2)^2+...
                               3/2*(x2*x1/(x2^2+y0^2)-x1^2/(x1^2+y0^2)+1/y0*x1*(atan(x2/y0)-atan(x1/y0))));
                        yConst5=((x2-x1)*y0/(x1^2+y0^2)-(atan(x2/y0)-atan(x1/y0)))/2;     
                        yConst10=((x1-x2)*y0/(x2^2+y0^2)+(atan(x2/y0)-atan(x1/y0)))/2;
                        y2Const4=x2/2*((x2/(x2^2+y0^2)-x1/(x1^2+y0^2))+(atan(x2/y0)-atan(x1/y0))/y0)...
                                 +(y0^2/(x2^2+y0^2)-y0^2/(x1^2+y0^2))/2; 
                        y2Const9=-x1/2*((x2/(x2^2+y0^2)-x1/(x1^2+y0^2))+(atan(x2/y0)-atan(x1/y0))/y0)...
                                  -(y0^2/(x2^2+y0^2)-y0^2/(x1^2+y0^2))/2;   
                        limConst2=x2/y0*(atan(x2/y0)-atan(x1/y0))-(log(x2^2+y0^2)-log(x1^2+y0^2))/2;
                        limConst7=-x1/y0*(atan(x2/y0)-atan(x1/y0))+(log(x2^2+y0^2)-log(x1^2+y0^2))/2;
                        limConst3=x2/2*(1/y0*(atan(x2/y0)-atan(x1/y0))-(x2/(x2^2+y0^2)-x1/(x1^2+y0^2)))...
                                  -(log((x2^2+y0^2)/(x1^2+y0^2))+y0^2*(1/(x2^2+y0^2)-1/(x1^2+y0^2)))/2;
                        limConst8=-x1/2*(1/y0*(atan(x2/y0)-atan(x1/y0))-(x2/(x2^2+y0^2)-x1/(x1^2+y0^2)))...
                                 +(log((x2^2+y0^2)/(x1^2+y0^2))+y0^2*(1/(x2^2+y0^2)-1/(x1^2+y0^2)))/2;

                       ywer1=a*yConst5+b*y2Const4;
                       ywer2=a*yConst10+b*y2Const9;
                       ywer3=c*yConst5+a*y2Const4;
                       ywer4=c*yConst10+a*y2Const9;
                       ymj1=a^3*(x2ycst1-ycst4)+3*a^2*b*(x2y2cst2-y2cst1)+3*a*b^2*(x2y3cst3-y3cst2)+b^3*(x2y4h2-y4cst3);
                       ymj2=a^3*(ycst4-x1ycst1)+3*a^2*b*(y2cst1-x1y2cst2)+3*a*b^2*(y3cst2-x1y3cst3)+b^3*(y4cst3-x1y4h2);
                       ymj3=a^2*c*(x2ycst1-ycst4)+(a^3+2*a*b*c)*(x2y2cst2-y2cst1)+(b^2*c+2*a^2*b)*(x2y3cst3-y3cst2)+b^2*a*(x2y4h2-y4cst3);
                       ymj4=a^2*c*(ycst4-x1ycst1)+(a^3+2*a*b*c)*(y2cst1-x1y2cst2)+(b^2*c+2*a^2*b)*(y3cst2-x1y3cst3)+b^2*a*(y4cst3-x1y4h2);
                       ymj5=a*c^2*(x2ycst1-ycst4)+(2*a^2*c+b*c^2)*(x2y2cst2-y2cst1)+(a^3+2*a*b*c)*(x2y3cst3-y3cst2)+a^2*b*(x2y4h2-y4cst3);
                       ymj6=a*c^2*(ycst4-x1ycst1)+(2*a^2*c+b*c^2)*(y2cst1-x1y2cst2)+(a^3+2*a*b*c)*(y3cst2-x1y3cst3)+a^2*b*(y4cst3-x1y4h2);
                       ymj7=c^3*(x2ycst1-ycst4)+3*c^2*a*(x2y2cst2-y2cst1)+3*c*a^2*(x2y3cst3-y3cst2)+a^3*(x2y4h2-y4cst3);
                       ymj8=c^3*(ycst4-x1ycst1)+3*c^2*a*(y2cst1-x1y2cst2)+3*c*a^2*(y3cst2-x1y3cst3)+a^3*(y4cst3-x1y4h2);
                       yme1=ywer1-4*ymj1;
                       me2=a^2*limConst3+2*a*b*yConst5+b^2*y2Const4;
                       yme3=ywer2-4*ymj2;
                       me4=a^2*limConst8+2*a*b*yConst10+b^2*y2Const9;
                       yme5=(1-2*NU)*ywer3-4*ymj3;
                       me6=a*c*limConst3+(a^2+b*c)*yConst5+a*b*y2Const4;
                       yme7=(1-2*NU)*ywer4-4*ymj4;
                       me8=a*c*limConst8+(a^2+b*c)*yConst10+a*b*y2Const9;
                       yme9=NU*ywer3-4*ymj3;
                       yme10=NU*ywer4-4*ymj4;
                       yme11=NU*ywer1-4*ymj5;
                       me12=c^2*limConst3+2*a*c*yConst5+a^2*y2Const4;
                       yme13=NU*ywer2-4*ymj6;
                       me14=c^2*limConst8+2*a*c*yConst10+a^2*y2Const9;
                       yme15=(1-2*NU)*ywer1-4*ymj5;
                       yme16=(1-2*NU)*ywer2-4*ymj6;
                       yme17=ywer3-4*ymj7;
                       yme18=ywer4-4*ymj8;
                       S111q1=2*yme1+2*n1*me2+n1*limConst2;
                       S111q2=2*yme3+2*n1*me4+n1*limConst7;
                       S211q1=2*yme5+4*NU*n1*me6+(1-2*NU)*2*n2*me2-(1-4*NU)*n2*limConst2;
                       S211q2=2*yme7+4*NU*n1*me8+(1-2*NU)*2*n2*me4-(1-4*NU)*n2*limConst7;
                       S112q1=2*yme9+(2-2*NU)*n1*me6+2*NU*n2*me2+(1-2*NU)*n2*limConst2;
                       S112q2=2*yme10+(2-2*NU)*n1*me8+2*NU*n2*me4+(1-2*NU)*n2*limConst7;
                       S212q1=2*yme11+2*n1*NU*me12+(2-2*NU)*n2*me6+(1-2*NU)*n1*limConst2;
                       S212q2=2*yme13+2*n1*NU*me14+(2-2*NU)*n2*me8+(1-2*NU)*n1*limConst7;
                       S122q1=2*yme15+4*NU*n2*me6+(1-2*NU)*2*n1*me12-(1-4*NU)*n1*limConst2;
                       S122q2=2*yme16+4*NU*n2*me8+(1-2*NU)*2*n1*me14-(1-4*NU)*n1*limConst7;
                       S222q1=2*yme17+2*n2*me12+n2*limConst2;
                       S222q2=2*yme18+2*n2*me14+n2*limConst7;
                       AS=Consts*[S111q1 S211q1 S111q2 S211q2;S112q1 S212q1 S112q2 S212q2;S122q1 S222q1 S122q2 S222q2]; 
                    end 
                    %和位移相关的体积力积分
                    GravC=-LE(ElNo)/(x2-x1)/8/pi/SG;
                    GravC1=((x2^2+y0^2)*log(x2^2+y0^2)-(x1^2+y0^2)*log(x1^2+y0^2))/2;
                    GravC2=y0*(x2*log(x2^2+y0^2)-x1*log(x1^2+y0^2)-2*((x2-x1)-y0*(atan(x2/y0)-atan(x1/y0)))+(x2-x1));
                    VForce=Den*AG*(n1*(a*GravC1+b*GravC2)+n2*(c*GravC1+a*GravC2))...
                          -1/2/(1-NU)*Den*AG'*[(a*GravC1+b*GravC2) (c*GravC1+a*GravC2)]'*[n1 n2]';
                    FV=VForce*GravC;
%                     体积力
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
%                 ----------------
                   AUi(:,1)=AUi(:,1)+Gi*[UPknow(En1,6:7) UPknow(En2,4:5)]';
                   APi(:,1)=APi(:,1)+Hi*[UPknow(En1,2:3) UPknow(En2,2:3)]';
                   ADi(:,1)= ADi(:,1)+AD*[UPknow(En1,6:7) UPknow(En2,4:5)]';
                   ASi(:,1)= ASi(:,1)+AS*[UPknow(En1,2:3) UPknow(En2,2:3)]';
                   FVjj(:,1)=FVjj(:,1)+FV;
                   SVjj(:,1)=SVjj(:,1)+[SV11 SV12 SV22]'; 
            end 
                iNterUi=(AUi-APi+FVjj)';
                iNterPi=(ADi-ASi+SVjj)';
                iNterU=cat(1,iNterU,[iA InNode(iKP,:) iNterUi]);
                iNterP=cat(1,iNterP,[iA InNode(iKP,:) iNterPi]);
               break
    end
        
end


if InCount~=0
    disp('internal point:');
    disp(sprintf('\t%s\t\t%s\t\t\t\t%s\t\t\t\t\t\t%s\t\t\t\t\t\t\t%s','Area','CX','CY','UX','UY'));
    disp(sprintf('%5.0f\t%12.8f\t%12.8f\t%24.12f\t%24.12f\n',iNterU'));
    disp(sprintf('\t%s\t\t%s\t\t\t\t%s\t\t\t\t\t%s\t\t\t\t\t%s\t\t\t\t\t%s','Area','CX','CY','SX','SXY','SY'))
    disp(sprintf('%5.0f\t%12.8f\t%12.8f\t%16.10f\t%16.10f\t%16.10f\n',iNterP'))
end

toc