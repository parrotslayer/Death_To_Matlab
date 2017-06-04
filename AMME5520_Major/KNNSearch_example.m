% KNNSearch example
clear all
clc

load fisheriris
x = meas(:,3:4);
y = [5 1.45;6 2;2.75 .75];

[n,d]=knnsearch(x,y,'k',10,'distance','minkowski','p',5);
[ncb,dcb] = knnsearch(x,y,'k',10,...
   'distance','chebychev');

gscatter(x(:,1),x(:,2),species)
line(y(:,1),y(:,2),'marker','x','color','k',...
   'markersize',10,'linewidth',2,'linestyle','none')
line(x(n,1),x(n,2),'color',[.5 .5 .5],'marker','o',...
   'linestyle','none','markersize',10)
line(x(ncb,1),x(ncb,2),'color',[.5 .5 .5],'marker','p',...
   'linestyle','none','markersize',10)
legend('setosa','versicolor','virginica','query point',...
'minkowski','chebychev','Location','best')