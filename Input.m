function  [ProblemType,Mat,Node,Area,Element,BCU,InNode]=INput(filename)
 % N  NUMBERS OF BOUNDARY ELEMENTS
 % L INTERNAL point NUMBERS(needed calculated)
 % E elastic module
 % XNU possion ratio
 
 fid=fopen(filename,'r');
 %Problem type? 1:plane strain,2:plane stress？
 line=fgetl(fid);
 ProblemType=fscanf(fid,'%g',[1 1]);
 
 %Material number:density, elastic modulus, Possion's ratio
 line=fgetl(fid); line=fgetl(fid);
 MatNo=fscanf(fid,'%g',[1 1]);
 Mat=fscanf(fid,'%g',[3 MatNo])';
 Mat(:,4)=Mat(:,2)/2./(1+Mat(:,3));
 
 %Node Coordinate
 line=fgetl(fid); line=fgetl(fid);
 NodeNo=fscanf(fid,'%g',[1 1]);
 Node=fscanf(fid,'%g',[3 NodeNo])';
 
%  Area no, material no of area
line=fgetl(fid); line=fgetl(fid);
 Area=fscanf(fid,'%g',[2 inf])';
 
 %Area No,Element No,Node1 No,Node2 No
 line=fgetl(fid);
 Element=fscanf(fid,'%g',[4 inf])';
 [B,IX]=sort(Element(:,2));          %暂不考虑多域问题
 Element=Element(IX,:);
 %Displacement of Node:Node No, type, value
 %types are: 1 = x, 2 = y, 3 = x = y
 line=fgetl(fid); 
 line=fgetl(fid); 
 BCU=fscanf(fid,'%g',[3 inf])';
 %Node in area to be calculated
 line=fgetl(fid);
 InNode=fscanf(fid,'%g',[2 inf])';

 status=fclose(fid);
 
 fid=fopen(strcat(filename(1:findstr(filename,'.')),'lis'),'r');
 ElNoPreALL=zeros(0,8);
 while ~feof(fid)
    while isempty(strfind(line,'LKEY')) & ~feof(fid)
        line=fgetl(fid);
    end
    ElNoPre=fscanf(fid,'%g',[8 inf])';
    line='';
    ElNoPreALL=cat(1,ElNoPreALL,ElNoPre);
 end 
 status=fclose(fid);
 
 ElPnPt=zeros(size(ElNoPreALL,1),5);
 ElPnPt(:,1)=ElNoPreALL(:,1);
 ElPnPt(find(ElNoPreALL(:,2)==1),[2 4])=ElNoPreALL(find(ElNoPreALL(:,2)==1),[3 7]);
 ElPnPt(find(ElNoPreALL(:,2)==2),[3 5])=ElNoPreALL(find(ElNoPreALL(:,2)==2),[3 7]);
 
 DL=Node(Element(:,4),2:end)-Node(Element(:,3),2:end);
 LE=sum(DL.*DL,2).^0.5;
 DL=DL./[LE LE];
 NE=DL*[0 -1;1 0];

 NXN_Y=[NE(:,1) -NE(:,2)];
 NYNX=[NE(:,2) NE(:,1)];
 N_XNY=[-NE(:,1) NE(:,2)];
 Element(:,5:8) = 1;
 Element(:,9:12) = 0;
 for iE=1:size(ElPnPt,1)
     iENo=ElPnPt(iE,1);
     iER=find(Element(:,2)==iENo);     
     Element(iER,9:12) =[sum(ElPnPt(iE,2:3).*NXN_Y(iER,:),2) sum(ElPnPt(iE,2:3).*NYNX(iER,:),2)...
     sum(ElPnPt(iE,4:5).*NXN_Y(iER,:),2) sum(ElPnPt(iE,4:5).*NYNX(iER,:),2)];
 end 

 for IE1=1:size(Element,1)
     NodeNo1=Element(IE1,3);
     NodeNo2=Element(IE1,4);
     Temp=BCU(find(BCU(:,1)==NodeNo1),:);
     UxNodeNo1=Temp(find(Temp(:,2)==1),1);
     UyNodeNo1=Temp(find(Temp(:,2)==2),1);
     Temp=BCU(find(BCU(:,1)==NodeNo2),:);
     UxNodeNo2=Temp(find(Temp(:,2)==1),1);
     UyNodeNo2=Temp(find(Temp(:,2)==2),1);
     if ~isempty(UxNodeNo1) & ~isempty(UyNodeNo1) & ~isempty(UxNodeNo2) & ~isempty(UyNodeNo2)
         Element(IE1,5:8)=0;     %单元边界全约束
     end
     if  ~isempty(UxNodeNo1) & ~isempty(UxNodeNo2)  %& abs(NE(IE1,:)*[0 1]')<1e-6
         Element(IE1,[5 7])=0;      %单元边界法向（X向）约束
     end
     if  ~isempty(UyNodeNo1) & ~isempty(UyNodeNo2) %& abs(NE(IE1,:)*[1 0]')<1e-6 
         Element(IE1,[6 8])=0;     %单元边界法向（Y向）约束
     end     
 end


 