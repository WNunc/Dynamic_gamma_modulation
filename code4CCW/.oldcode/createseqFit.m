function [fitresult, gof] = createseqFit(ydata, xdata, wdata)
%CREATEFIT(YDATA,XDATA,WDATA)
%  Create a thetaseqence fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : ydata ==>tbins of theta sequence
%      Y Output: xdata ==>xbins
%      Weights : wdata ==>pxn closest pxn at COM or peak probability
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 07-Mar-2022 22:38:32 自动生成


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( ydata, xdata);

% Set up fittype and options.
ft = fittype( 'poly1' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'xdata vs. ydata with wdata', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'ydata', 'Interpreter', 'none' );
% ylabel( 'xdata', 'Interpreter', 'none' );
% grid on


