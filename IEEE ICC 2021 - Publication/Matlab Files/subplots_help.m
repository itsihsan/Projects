clear;clc;clf
fig = figure;

n = -10:0.01:10;

subplot(3,1,1);
a = logsig(n);
plot(n,a,'b','LineWidth',1.5)
ylim([-0.1 1.1])
legend('$y = \frac{1}{1+e^{-x}}$','FontSize',12,'location','northwest','Interpreter','latex')
%legend('Sigmoid activation function','location','northwest')
ylabel('Output (y)');
xlabel('Input (x)');
hold on

subplot(3,1,2);
b = tansig(n);
plot(n,b,'r','LineWidth',1.5)
ylim([-1.1 1.1])
legend('$y = \frac{e^x-e^{-x}}{e^x+e^{-x}}$','FontSize',12,'location','northwest','Interpreter','latex')
%legend('Tangent sigmoid activation function','location','northwest')
ylabel('Output (y)');
xlabel('Input (x)');

subplot(3,1,3);
c = poslin(n);
plot(n,c,'k','LineWidth',1.5)
ylim([-1 10])
legend('$y = max(0,x)$','FontSize',12,'location','northwest','Interpreter','latex')
ylabel('Output (y)');
xlabel('Input (x)');
%legend('Positive linear activation function','location','northwest')

% Give common xlabel, ylabel and title to your figure
han=axes(fig,'visible','off'); 
% han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Activation functions output value');
xlabel(han,'Activation functions input value');
% title(han,'yourTitle');


% Below help is for bar graph\
% y = [15,28,33;30,52,56;62,84,85]; %5N-Half
% y = [30,31,38;58,57,63;89,87,90]; %5N-Quarter
% y = [ 34,33;52,56;78,85]; %FP, 5N-8neurons-Half resolution
% y = [34,38;52,63;78,90]; %FP, 5N-8neurons-quarter resolution
% y = [33,17,34;58,41,59;88,81,87]; %4N-Half resolution
% y = [21,19;37,35;62,65]; 
y = [28.5,34.3;47.1,58.1;74.4,83.6]; 
x=[39,66,95]
% x=['1st','2nd','3rd']
b = bar(x,y)
b(1).FaceColor = [0.5 0.5 0.5];
% b(1).FaceColor = [0 0 1];
b(2).FaceColor = [0 1 0];
% b(2).FaceColor = [0.75 0.75 0.75];
% b(3).FaceColor = [0 1 0];

xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = string(b(2).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

xtips3 = b(3).XEndPoints;
ytips3 = b(3).YEndPoints;
labels3 = string(b(3).YData);
text(xtips3,ytips3,labels3,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

legend('RSS based algorithm','NNLEF $(P_n = 8)$','FontSize',11,'location','northwest','Interpreter','latex')
% legend('Fingerpinting framework','Framework-C','location','northwest')

% xlabel({'Confidence Ellipse Number'},'Interpreter','latex','FontSize',11)
xlabel({'Confidence Ellipses with Different CLs (%)'},'Interpreter','latex','FontSize',12)
ylabel({'\% of Estimated Locations'},'Interpreter','latex','FontSize',12)


grid on
grid minor

