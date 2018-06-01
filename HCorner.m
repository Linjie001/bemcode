function [H1 H2]=HCorner(Nu,Ang1,Ang2)
if Ang2<Ang1
    Ang2=Ang2+2*pi
end
Ang21=Ang2-Ang1

Ang2m=Ang1*2+Ang21   %Òª¼ì²é
c0=-1/(8*pi*(1-Nu)*sin(Ang21))
c1=(1-2*Nu)*(1-cos(Ang21))
c2=1-cos(Ang21)+0.5*(sin(2*Ang1)*sin(Ang21)-Ang21*sin(Ang2m))
c3=0.5*(Ang21*cos(Ang2m)-cos(2*Ang1)*sin(Ang21))
c4=(1-2*Nu)*(sin(Ang21)*(2*log(2*sin(Ang21/2))-1)-Ang21*cos(Ang21))
c5=c1
c6=1-cos(Ang21)+0.25*(cos(2*Ang2)-cos(2*Ang1)+2*Ang21*sin(2*Ang1))
c7=0.25*(sin(2*Ang2)-sin(2*Ang1)-2*Ang21*cos(2*Ang1))
c8=(1-2*Nu)*(Ang21+sin(Ang21))
c9=2*(1-cos(Ang21))
H1=c0*[c1+c2 c3+c4    c5+c6 c7+c8;
       c3-c4 c1-c2+c9 c7-c8 c5-c6+c9]
d0=c0
d1=c1
d2=c9-c6
d3=-c7
d4=-c8
d5=c1
d6=1-cos(Ang21)-0.5*(sin(2*Ang2)*sin(Ang21)-Ang21*sin(Ang2m))
d7=0.5*(cos(2*Ang2)*sin(Ang21)-Ang21*cos(Ang2m))
d8=-c4
d9=c9
H2=d0*[d1+d2 d3+d4    d5+d6 d7+d8;
       d3-d4 d1-d2+d9 d7-d8 d5-d6+d9] 

