%-------------------------------------------------------------------
% polyfit - order = ord
ord = 15;
% xscale = 0:50:400;
xscale = 1:1:575;
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

        if c_f(tt,ii) > 1
            c_f(tt,ii) = 1;
        end

        if c_f(tt,ii) < 0
            c_f(tt,ii) = 0;
        end
    end
end
ihsan = figure(2)
% saveas (figure,'Total Error','epsc')
axes1 = axes;
hold(axes1,'on');
xscale_adj = xscale;

g(1)=plot(xscale_adj,c_f(1,1:1:count),'r-o','MarkerIndices',90:90:length(c),'MarkerFaceColor','w','LineWidth',1.5);
hold on;
g(2)=plot(xscale_adj,c_f(2,1:1:count),'b-s','MarkerIndices',90:90:length(c),'MarkerFaceColor','w','LineWidth',1.5);
hold on;
g(3)=plot(xscale_adj,c_f(3,1:1:count),'k-d','MarkerIndices',90:90:length(c),'MarkerFaceColor','w','LineWidth',1.5);
hold on;
g(4)=plot(xscale_adj,c_f(4,1:1:count),'r-.o','MarkerIndices',20:50:length(c),'MarkerFaceColor','w','LineWidth',1.5);
hold on;
g(5)=plot(xscale_adj,c_f(5,1:1:count),'b-.s','MarkerIndices',30:50:length(c),'MarkerFaceColor','w','LineWidth',1.5);
hold on;
g(6)=plot(xscale_adj,c_f(6,1:1:count),'k-.d','MarkerIndices',40:50:length(c),'MarkerFaceColor','w','LineWidth',1.5);
hold on;

xlabel('^{d}/_{d}')
xlabel('d_r/d_i','FontSize',11)
xlabel('$\frac{d_r}{d_i}$','Interpreter','latex','FontSize',19)
xlabel('FontSize',11)
ylabel('RSS_i (dBm)','FontSize',11)
legend('r = 50m','r = 75m','r = 100m','FontSize',11)

% Rounding off y-axis to 2 decimal points
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))

grid on
set(gca,'GridColor','[0 0 0]')
set(gca,'GridLineStyle','-.')
% set(gca,'GridColor','[0 0 0]')
xlim([0 3000])
ylim([0.08 0.3])

% set(gca,'Color','[1 1 1]') % setting the background color of the plot

set(gca,'Color','[0.88 0.88 0.88]') % setting the background color of the plot
set(gca,'GridColor','[0.5 0.5 0.5]')
% set(gca,'GridColor','[0 0 0]')
set(gca, ...
'GridColor'   ,'[0 0 0]'  , ...
'Box'         , 'off'     , ...
'TickDir'     , 'out'     , ...
'TickLength'  , [.02 .02] , ...
'XMinorTick'  , 'on'      , ...
'YMinorTick'  , 'on'      , ...
'YGrid'       , 'on'      , ...
'XColor'      , [.3 .3 .3], ...
'XTick'       , 0:90:540, ...
'YColor'      , [.3 .3 .3], ...
'YTick'       , 0:0.1:0.6, ...
'LineWidth'   , 1         );

% 'XColor'      , [.3 .3 .3], ...
% 'XTick'       , 1.0:0.2:2.2, ...
% 'YColor'      , [.3 .3 .3], ...
% 'YTick'       , -85:5:-55, ...






xL = get(gca, 'XLim');
plot(xL, [0.18 0.18], '-r', 'LineWidth', 2.0)

% xL = get(gca, 'XLim');
% plot(xL, [0.16 0.16], '-r', 'LineWidth', 2.0)


xlabel('Training Data ','FontSize',11);
ylabel('Total Error','FontSize',11);
legend({'Po = 0.5','Po = 0.3','Po = 0.1','Po = 0.01','Po = 0.001','Po = 0.0005','Total Error - LRT'},'FontSize',11)
% legend({'NLoS 300ns','NLoS 500ns','NLoS 700ns'},'FontSize',11)

grid on
set(gca,'GridColor','[0.5 0.5 0.5]')
set(gca,'GridLineStyle','-.')

% Arrow Code

br = annotation('arrow');
d = br.Color;
br.Color = 'r';
set(br,'LineStyle','-.','LineWidth',1.5,'X',[0.43 0.13],'Y',[0.505 0.505])%non optimal
% set(br,'LineStyle','-.','LineWidth',1.5,'X',[0.43 0.13],'Y',[0.63 0.63])
cc=text(70,0.29,'LRT Total Error - 100m','Color','red','FontSize',11)


cr = annotation('arrow');
e = cr.Color;
cr.Color = 'b';
set(cr,'LineStyle','-.','LineWidth',1.5,'X',[0.43 0.13],'Y',[0.41 0.41]) % non optimal
% set(cr,'LineStyle','-.','LineWidth',1.5,'X',[0.43 0.13],'Y',[0.605 0.605])
bb=text(70,0.13,'LRT Total Error - 75m','Color','blue','FontSize',11)

ar = annotation('arrow');
c = ar.Color;
ar.Color = 'k';
set(ar,'LineStyle','-.','LineWidth',1.5,'X',[0.43 0.13],'Y',[0.18 0.18]) % non optimal
% set(ar,'LineStyle','-.','LineWidth',1.5,'X',[0.43 0.13],'Y',[0.50 0.50])
aa=text(70,0.072,'LRT Total Error - 50m','Color','black','FontSize',11)

% xlim([0 4000])
% legend('NLoS  50 (ns)','NLoS 100 (ns)','NLoS 150 (ns)','NLoS 200 (ns)','NLoS 250 (ns)','NLoS 300 (ns)','NLoS 350 (ns)','TLB','FontSize',11)
% legend('Thermal Noise 300 (ns)','Thermal Noise 200 (ns)','Thermal Noise 100 (ns)','Thermal Noise  50  (ns)','Thermal Noise  10  (ns)','Thermal Noise  1    (ns)','FontSize',11)
% legend('NLoS  50  (ns)','NLoS 100 (ns)','NLoS 150 (ns)','NLoS 200 (ns)','NLoS 250 (ns)','NLoS 300 (ns)','NLoS 350 (ns)','TLB','FontSize',11)
% grid on
% set(gca,'GridColor','[0.5 0.5 0.5]')
% set(gca,'GridLineStyle','-.')
% xlabel('Deployement Data - NLoS (ns) ','FontSize',11);
% legend('NN-LT NLoS  50 (ns)','NN-LT NLoS 100 (ns)','NN-LT NLoS 150 (ns)','NN-LT NLoS 200 (ns)','NN-LT NLoS 250 (ns)','FontSize',11)
% legend('NN-LT NLoS  50  (ns)','NN-LT NLoS 100 (ns)','NN-LT NLoS 150 (ns)','NN-LT NLoS 200 (ns)','NN-LT NLoS 250 (ns)','FontSize',11)
% legend('NN-LT NLoS  50ns','NN-LT NLoS 100ns','NN-LT NLoS 150ns','NN-LT NLoS 200ns','NN-LT NLoS 250ns','FontSize',11)
% legend('NN-LT NLoS  50 ns','NN-LT NLoS 100 ns','NN-LT NLoS 150 ns','NN-LT NLoS 200 ns','NN-LT NLoS 250 ns','FontSize',11)
% legend('NLoS  50  ns','NLoS 100 ns','NLoS 150 ns','NLoS 200 ns','NLoS 250 ns','NLoS 300 ns','NLoS 350 ns','TLB','FontSize',11)
% legend('Thermal Noise 300 ns','Thermal Noise 200 ns','Thermal Noise 100 ns','Thermal Noise  50  ns','Thermal Noise  10  ns','Thermal Noise  1    ns','FontSize',11)
% legend('NN-LT NLoS  50ns','NN-LT NLoS 100ns','NN-LT NLoS 150ns','NN-LT NLoS 200ns','NN-LT NLoS 250ns','FontSize',11)
% legend('NLoS  50ns','NLoS 100ns','NLoS 150ns','NLoS 200ns','NLoS 250ns','NLoS 300ns','NLoS 350ns','TLB','FontSize',11)
% legend('Thermal Noise 300ns','Thermal Noise 200ns','Thermal Noise 100ns','Thermal Noise  50ns','Thermal Noise  10ns','Thermal Noise  1ns','FontSize',11)
% legend('r = 200 meters','r = 175 meters','r = 150 meters','r = 125 meters','FontSize',11)