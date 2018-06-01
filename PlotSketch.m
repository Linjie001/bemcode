function PlotSketch(ScalarFactor,NodeU)
global Element  Node
hold on
XYMin=min(Node(:,2:3));
XYMax=max(Node(:,2:3));
DXYMax=XYMax-XYMin;
XYMin=XYMin-DXYMax/20;
XYMax=XYMax+DXYMax/20;
axis([XYMin(1) XYMax(1) XYMin(2) XYMax(2)])
axis equal
for iE=1:size(Element,1)
    Nab=Element(iE,3:4);
    x=Node(Nab,2)+ScalarFactor*NodeU(Nab,2);
    y=Node(Nab,3)+ScalarFactor*NodeU(Nab,3);
    if ScalarFactor==0
        plot(x,y);
    else
        plot(x,y,'r');
    end
end
hold off
