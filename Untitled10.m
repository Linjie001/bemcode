clc
clear
format long

q=10;
% A=[0.707107   0.7072
% 0.707108   0.7072
% 0.70711    0.7072
% 0.707115   0.7072
% 0.70712   0.7072];
% B=[   9.5619051298	   -7.9945038529	    9.5637432874
%   9.5618990083	   -7.9944977217	    9.5637492728
%   9.5618867650	   -7.9944854593	    9.5637612436
%   9.5618561560	   -7.9944548034	    9.5637911717
%   9.5618255457	   -7.9944241476	    9.5638211009
%  10.2578640311	   -0.5318441151	   10.2564896943
% ];
% C=[157.1696	        157.4200	        157.9941 
%  157.1645	        157.4249	        157.9804 
%  157.1542	        157.4348	        157.9528 
%  157.1287	        157.4596	        157.8838 
%  157.1031	        157.4844	        157.8148 ];
A=[70.710         70.710
70.710005      70.710
70.710         70.710
70.710010      70.710
70.71005         70.710
70.71010     70.710
];
B=[ 10.0002933560	    0.0005439060	   10.0003777984 
 10.0002933041	    0.0005439607	   10.0003785732 
 10.0002933560	    0.0005439060	   10.0003777984 
 10.0002932060	    0.0005440158	   10.0003793973 
 10.0002926198	    0.0005443277	   10.0003857625 
 10.0002915025	    0.0005444048	   10.0003940533 
];
C=[   -2042.0374	         -4.9075	      -2042.0374
  -2042.0223	         -4.9226	      -2042.0376
  -2042.0374	         -4.9075	      -2042.0374
  -2042.0072	         -4.9377	      -2042.0379
  -2041.8857	         -5.0592	      -2042.0403
  -2041.7325	         -5.2125	      -2042.0449
];
r=sqrt(A(:,1).^2+A(:,2).^2);
u=atan(A(:,2)./A(:,1));
m=[];
n=[];
p=[];
q=[];
for ih=1:size(u,1)
    
%     a=[(cos(u(ih)))^2 (sin(u(ih)))^2 -sin(2*u(ih));...
%         (sin(u(ih)))^2 (cos(u(ih)))^2 sin(2*u(ih));...
%         sin(2*u(ih))/2  -sin(2*u(ih))/2  cos(2*u(ih))]*[sig1 sig2 0]';
    a=(cos(u(ih)))^2 *B(ih,1)+(sin(u(ih)))^2*B(ih,3)+sin(2*u(ih))*B(ih,2);
    b=(sin(u(ih)))^2 *C(ih,1)+(cos(u(ih)))^2*C(ih,3)-sin(2*u(ih))*C(ih,2);
    sig1=10*(1-1/r(ih)^2);
%     sig2=10*(1+1/r(ih)^2);
              m=cat(2,m,a);
              q=cat(2,q,b);
              n=cat(2,n,sig1);
%               p=cat(2,p,sig2);
end
m'
n'
 q'
% p
