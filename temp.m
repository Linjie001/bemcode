clc
clear
A=[ 2 2 
    ];
x=A(:,1);
y=A(:,2);
dis=sqrt(x.^2+y.^2);
u=atan(y./x);
aq=[];
for ij=1:size(u,1)
    trans=[cos(u(ij)) -sin(u(ij));sin(u(ij)) cos(u(ij))];
    sig1=10*(1-dis^2/1);
    sig2=10*(1+dis^2/1);
    sig=trans*[sig1 sig2]';
    aq=cat(1,aq,sig(1));
end
aq

