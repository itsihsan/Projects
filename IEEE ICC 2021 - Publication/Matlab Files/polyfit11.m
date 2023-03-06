clf
%-------------------------------------------------------------------
% polyfit - order = ord
ord = 7;
% xscale = 0:50:400;
% xscale = 1:1:400;
count=length(xscale);
for i = 1:size(c,1)
    coefs(i,:) = polyfit(xscale,c(i,:),ord);
    c_f(i,:) = polyval(coefs(i,:),xscale);
end

c_f_original = c_f;
[row,col]=size(c)
%-------------------------------------------------------------------
for tt = 1:row
    for ii = 1:1:count

        if c_f(tt,ii) > 100
            c_f(tt,ii) = 100;
        end

        if c_f(tt,ii) < 0
            c_f(tt,ii) = 0;
        end
    end
end
ihsan = figure(11)
% saveas (figure,'Total Error','epsc')
axes1 = axes;
hold(axes1,'on');
xscale_adj = xscale;

g(1)=plot(xscale_adj,c_f(1,1:1:count),'g-','LineWidth',1.5);
hold on;
g(2)=plot(xscale_adj,c_f(2,1:1:count),'b-','LineWidth',1.5);
hold on;
g(3)=plot(xscale_adj,c_f(3,1:1:count),'k-','LineWidth',1.5);

sigma1 = get(gca, 'XLim');
plot(sigma1, [39 39], '-.g', 'LineWidth', 1.5)
sigma2 = get(gca, 'XLim');
plot(sigma2, [66 66], '-.b', 'LineWidth', 1.5)
sigma3 = get(gca, 'XLim');
plot(sigma3, [95 95], '-.k', 'LineWidth', 1.5)

xlabel({'Number of Neurons ($P_n$)'},'Interpreter','latex','FontSize',11)
% ylabel({'% of Estimated Locations'},'Interpreter','latex','FontSize',11)
ylabel({'\% of NNLEF Estimated Locations'},'Interpreter','latex','FontSize',11)
hlegend = legend('Est. locations - CL 39%','Est. locations - CL 66%','Est. locations - CL 95%','Bound - CL 39%','Bound - CL 66%','Bound - CL 95%','FontSize',11,'location','southeast','Interpreter','latex')
hlegend.NumColumns=2
% set(hlegend,'color','[0.95 0.95 0.95]');

set(gca,'GridColor','[0 0 0]')
set(gca,'GridLineStyle','-')
xlim([0 50])
ylim([0 100])
% set(gca,'Color','[0.95 0.95 0.95]')
set(gca,'Color','[1 1 1]')
set(gca, ...
'GridColor'   ,'[0 0 0]'  , ...
'Box'         , 'off'     , ...
'TickDir'     , 'out'     , ...
'TickLength'  , [.02 .02] , ...
'XMinorTick'  , 'on'      , ...
'YMinorTick'  , 'on'      , ...
'YGrid'       , 'on'      , ...
'XColor'      , [0 0 0], ...
'YColor'      , [0 0 0], ...
'XTick'       , 0:5:50, ...
'YTick'       , 0:10:100, ...
'LineWidth'   , 1         );

