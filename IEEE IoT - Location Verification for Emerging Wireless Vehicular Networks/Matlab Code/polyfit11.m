%-------------------------------------------------------------------
clf;
% polyfit - order = ord
ord = 7;
% xscale = 0:50:400;
xscale = 1:1:5210;
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
ihsan = figure(16)
% saveas (figure,'Total Error','epsc')
axes1 = axes;
hold(axes1,'on');
xscale_adj = xscale;

g(1)=plot(xscale_adj,c_f(1,1:1:count),'b-','MarkerIndices',100:700:length(c),'MarkerFaceColor','w','LineWidth',2);
hold on;
g(2)=plot(xscale_adj,c_f(2,1:1:count),'b--','MarkerIndices',100:700:length(c),'MarkerFaceColor','w','LineWidth',2);
hold on;
g(3)=plot(xscale_adj,c_f(3,1:1:count),'b-.','MarkerIndices',100:700:length(c),'MarkerFaceColor','w','LineWidth',2);
hold on;
g(4)=plot(xscale_adj,c_f(4,1:1:count),'b:','MarkerIndices',100:700:length(c),'MarkerFaceColor','w','LineWidth',2);
hold on;

% Arrow Code
ar = annotation('arrow');
b = ar.Color;
ar.Color = 'b';
set(ar,'LineStyle','-','LineWidth',2,'X',[0.43 0.13],'Y',[0.54 0.54])

br = annotation('arrow');
c = ar.Color;
br.Color = 'b';
set(br,'LineStyle','--','LineWidth',2,'X',[0.43 0.13],'Y',[0.445 0.445])

cr = annotation('arrow');
d = ar.Color;
cr.Color = 'b';
set(cr,'LineStyle','-.','LineWidth',2,'X',[0.43 0.13],'Y',[0.32 0.32])

dr = annotation('arrow');
e = ar.Color;
dr.Color = 'b';
set(dr,'LineStyle',':','LineWidth',2,'X',[0.43 0.13],'Y',[0.25 0.25])





g(1)=plot(xscale_adj,c_f(1,1:1:count),'k-','MarkerIndices',100:700:length(c),'MarkerFaceColor','w','LineWidth',2);
hold on;
g(2)=plot(xscale_adj,c_f(2,1:1:count),'r-','MarkerIndices',100:700:length(c),'MarkerFaceColor','w','LineWidth',2);
hold on;
g(3)=plot(xscale_adj,c_f(3,1:1:count),'g-','MarkerIndices',100:700:length(c),'MarkerFaceColor','w','LineWidth',2);
hold on;
g(4)=plot(xscale_adj,c_f(4,1:1:count),'m-','MarkerIndices',100:700:length(c),'MarkerFaceColor','w','LineWidth',2);
hold on;
g(5)=plot(xscale_adj,c_f(5,1:1:count),'b-','MarkerIndices',100:700:length(c),'MarkerFaceColor','w','LineWidth',2);
hold on;
xlim([0 5000])

g(6)=plot(xscale_adj,c_f(6,1:1:count),'k-.d','MarkerIndices',100:700:length(c),'MarkerFaceColor','w','LineWidth',1.5);
hold on;


legend('$C_{10}=0.4,C_{01}=0.6$','$C_{10}=0.3,C_{01}=0.7$','$C_{10}=0.2,C_{01}=0.8$','$C_{10}=0.1,C_{01}=0.9$','Interpreter','latex')

xlabel('Time in Seconds ','FontSize',11);
ylabel('False Positive Rate','FontSize',11);
legend({'r = 150m','r = 200m','r = 300m','r = 1000m','r = 10000m'},'FontSize',11)

xlabel('Time in Seconds ','FontSize',11);
ylabel('Detection Rate','FontSize',11);
legend({'r = 150m','r = 200m','r = 300m','r = 1000m','r = 10000m'},'FontSize',11)

%same colour
g(1)=plot(xscale_adj,c_f(1,1:1:count),'k-','MarkerIndices',100:700:length(c),'MarkerFaceColor','w','LineWidth',2);
hold on;
xlim([0 5000])
g(2)=plot(xscale_adj,c_f(2,1:1:count),'k--','MarkerIndices',100:700:length(c),'MarkerFaceColor','w','LineWidth',2);
hold on;
g(3)=plot(xscale_adj,c_f(3,1:1:count),'k-.','MarkerIndices',100:700:length(c),'MarkerFaceColor','w','LineWidth',2);
hold on;
g(4)=plot(xscale_adj,c_f(4,1:1:count),'k:','MarkerIndices',100:700:length(c),'MarkerFaceColor','w','LineWidth',2);
hold on;

% Arrow Code
ar = annotation('arrow');
b = ar.Color;
ar.Color = 'k';
set(ar,'LineStyle','-.','LineWidth',2,'X',[0.43 0.13],'Y',[0.42 0.42])

br = annotation('arrow');
c = ar.Color;
br.Color = 'r';
set(br,'LineStyle','-.','LineWidth',2,'X',[0.43 0.13],'Y',[0.34 0.34])

cr = annotation('arrow');
d = ar.Color;
cr.Color = 'g';
set(cr,'LineStyle','-.','LineWidth',2,'X',[0.43 0.13],'Y',[0.25 0.25])

dr = annotation('arrow');
e = ar.Color;
dr.Color = 'm';
set(dr,'LineStyle','-.','LineWidth',2,'X',[0.43 0.13],'Y',[0.19 0.19])

er = annotation('arrow');
f = ar.Color;
er.Color = 'b';
set(er,'LineStyle','-.','LineWidth',2,'X',[0.43 0.13],'Y',[0.2 0.2])

xlabel('Time in Seconds','FontSize',11)
ylabel('Total Error','FontSize',11)
legend('Log Sigmoid','Tangent Sigmoid','ReLu','Log Sigmoid','Tangent Sigmoid','ReLu','FontSize',11)

% Rounding off y-axis to 2 decimal points
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'))

grid on
set(gca,'GridColor','[0 0 0]')
set(gca,'GridLineStyle','-.')
% set(gca,'GridColor','[0 0 0]')
xlim([0 3000])
ylim([0.08 0.3])

set(gca,'Color','[0.88 0.88 0.88]')
set(gca, ...
'GridColor'   ,'[1 1 1]'  , ...
'Box'         , 'off'     , ...
'TickDir'     , 'out'     , ...
'TickLength'  , [.02 .02] , ...
'XMinorTick'  , 'on'      , ...
'YMinorTick'  , 'on'      , ...
'YGrid'       , 'on'      , ...
'XColor'      , [.3 .3 .3], ...
'YColor'      , [.3 .3 .3], ...
'YTick'       , 0.1:0.0273:0.40, ...
'LineWidth'   , 1         );

xL = get(gca, 'XLim');
plot(xL, [0.18 0.18], '-r', 'LineWidth', 2.0)

% xL = get(gca, 'XLim');
% plot(xL, [0.16 0.16], '-r', 'LineWidth', 2.0)

xlabel('Time in Seconds ','FontSize',11);
ylabel('False Positive Rate','FontSize',11);
legend({'r = 150m','r = 200m','r = 300m','r = 1000m','r = 10000m'},'FontSize',11)



xlabel('Time in Seconds ','FontSize',11);
ylabel('Total Error','FontSize',11);
legend({'Po = 0.5','Po = 0.3','Po = 0.1','Po = 0.01','Po = 0.001','Po = 0.0005','Total Error - LRT'},'FontSize',11)
% legend({'NLoS 300ns','NLoS 500ns','NLoS 700ns'},'FontSize',11)

grid on
set(gca,'GridColor','[0.5 0.5 0.5]')
set(gca,'GridLineStyle','-.')

% Arrow Code
ar = annotation('arrow');
c = ar.Color;
ar.Color = 'k';
set(ar,'LineStyle','-.','LineWidth',2,'X',[0.43 0.13],'Y',[0.87 0.87])

br = annotation('arrow');
d = br.Color;
br.Color = 'r';
set(br,'LineStyle','-.','LineWidth',2,'X',[0.5 0.13],'Y',[0.72 0.72])

cr = annotation('arrow');
e = cr.Color;
cr.Color = 'b';
set(cr,'LineStyle','-.','LineWidth',2,'X',[0.5 0.13],'Y',[0.79 0.79])


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