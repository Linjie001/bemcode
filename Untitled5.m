clc
clear
format long
E=5;
v=0.3;
A=[4.7072    0.7072
2.7072    0.7072
1.1072    0.7072
0.9072    0.7072
0.7472    0.7072
0.7272    0.7072
0.7112    0.7072
0.7092    0.7072];
% s1代表数值解sx sy sxy
% s2代表边界元解析解
x=A(:,1);
y=A(:,2);
dis=sqrt(x.^2+y.^2);
u=atan(y./x);
aq=[];
for ij=1:size(u,1)
    trans=[cos(u(ij)) -sin(u(ij));sin(u(ij)) cos(u(ij))];
    UP=10*(1-v^2)/E*(1/(1-v)/dis(ij)+(1-2*v)/(1-v)*dis(ij));
    UQ=0;
    up=trans*[UP UQ]';
    aq=cat(1,aq,up(1));
end
aq