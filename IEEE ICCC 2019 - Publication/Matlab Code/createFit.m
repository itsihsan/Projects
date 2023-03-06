function [fitresult, gof] = createFit(x, y)
%CREATEFIT(X2,Y2)
%  Create a fit.
%
%  Data for 'Ihsan' fit:
%      X Input : x2
%      Y Output: y2
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 13-Mar-2019 16:45:53


%% Fit: 'Ihsan'.
[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( 'a*10*log10(x)+b', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.140691093663118 0.886794158592036];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'Ihsan3' );
h = plot( fitresult, xData, yData );
% legend('RSS Measurments', 'Curve Fitting', 'Location', 'NorthEast','FontSize',11 );
% legend(h,'RSS measurements', '$RSS_i=RSS_{d_r}+10\,\gamma\,log_{10}\Big(\frac{d_r}{d_i}\Big)$','Interpreter','Latex','Location', 'NorthWest','FontSize',11 );
% Label axes
%xlabel x
%ylabel y
% xlabel('$\frac{d_r}{d_i}$','Interpreter','latex','FontSize',19)
xlabel('Logrithmic Distance','FontSize',11)
ylabel('RSS','FontSize',11)

grid on

