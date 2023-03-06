y1 = smooth(y(1,:));
y2 = smooth(y(2,:));
y3 = smooth(y(3,:));
y4 = smooth(y(4,:));
y5 = smooth(y(5,:));
y6 = smooth(y(6,:));
y7 = smooth(y(7,:));
% y2 = y(2,:);
% y3 = y(3,:);
% y4 = y(4,:);
% y5 = y(5,:);
% y6 = y(6,:);
% plot(x,y1,'r-','LineWidth',1.5);
% hold on;
% plot(x,y2,'b-','LineWidth',1.5);
% plot(x,y3,'k-','LineWidth',1.5);
% plot(x,y4,'m-','LineWidth',1.5);
% plot(x,y5,'c-','LineWidth',1.5);
% plot(x,y6,'g-','LineWidth',1.5);
% plot(x,y7,'y-','LineWidth',1.5);
plot(x,y1,'Color','[1 0 0]','LineWidth',1.5);
hold on;
plot(x,y2,'Color','[0 0 1]','LineWidth',1.5);
plot(x,y3,'Color','[0 0 0]','LineWidth',1.5);
plot(x,y4,'Color','[1 0 1]','LineWidth',1.5);
plot(x,y5,'Color','[1 0.6471 0]','LineWidth',1.5);
plot(x,y6,'g-','LineWidth',1.5);
% plot(x,y7,'y-','LineWidth',1.5);
plot(x,y7,'Color',[0.6039 0.8039 0.1961],'LineWidth',1.5);

% plot(x,y)
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
'XTick'       , 0:20:120, ...
'YTick'       , 0:300:1800, ...
'LineWidth'   , 1         );
xticks([0 3 6 9 12 15 18])
% xticks([0 5 10 15 20])
yticks([0 4 8 12 16 20])
xticklabels({'0','3','6','9','12','15','18'})
yticklabels({'0','8','16','24','32','40'})

% xticklabels({'0 to 10','11 to 20','21 to 30','31 to 40','41 to 50','51 to 60','61 to 70','71 to 80','81 to 90','91 to 100','101 to 110','111 to 120','121 to 130','131 to 140','141 to 150','151 to 160','161 to 170','171 to 180','181 to 190','191 to 200'},'Interpreter','latex')
xticklabels({'0 to 10','11 to 20','21 to 30','31 to 40','41 to 50','51 to 60','61 to 70','71 to 80','81 to 90','91 to 100','101 to 110','111 to 120','121 to 130','131 to 140','141 to 150','151 to 160','161 to 170','171 to 180','181 to 190','191 to 200'})
xtickangle(90)
xlabel({'Mean Square Error (meters)'},'Interpreter','latex','FontSize',13)
ylabel({'Number of Samples'},'Interpreter','latex','FontSize',13)
hlegend = legend('$\,\ 3$ Neurons','$\,\ 8$ Neurons','20 Neurons','30 Neurons','40 Neurons','50 Neurons','Interpreter','latex','FontSize',12,'location','north')
hlegend.NumColumns=2
