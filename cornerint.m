syms st1 st2 st
st21=st2-st1
fc1=sin(st21-st)
c1=int(fc1,st,0,st21)
fc2=2*(sin(st1+st/2)).^2*sin(st21-st)
c2=int(fc2,st,0,st21)
fc3=-sin(2*st1+st)*sin(st21-st)
c3=int(fc3,st,0,st21)
fc4=sin(st21-st)/tan(st/2)
c4=int(fc4,st,0,st21)
pretty(c4)