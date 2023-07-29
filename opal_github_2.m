clear all; close all; clc; load opal.mat; % clear workspace & load relevant data
b = .65; k = 1.02; g = .92; % regional means of parameters other than alpha_Si
beta = -((O.*Z.^(-b)./k).^(1/g)-I)./S; % invert for ratio of alphas
basin = ones(size(X)); basin(Y<-30) = 1; basin(abs(Y)<30.001 & X>-80 & X<30) = 5; basin(Y>30 & X>-80 & X<30) = 4; basin(abs(Y)<30.001 & X<-80 | X>30 & abs(Y)<30.001) = 3; basin(Y>30 & X<-80 | X>30 & Y>30) = 2; % define regions for plot color
scatter(Si80,k.^g.*beta,10,basin,'filled') % plot
set(gca,'xscale','log','yscale','log','ticklabelinterpreter','latex','fontsize',16)
axis([.8 100 .0004 150])
caxis([0 6])
box on;
ylabel('$\alpha_{Si}$ (g OC/g Si)','interpreter','latex')
xlabel('[Si] at 80m ($\mu$mol/kg)','interpreter','latex')
hold on;
scatter([1 1 1 1 1],[.01 .005 .0025 .00125 .000625],50,[1 2 3 4 5],'filled')
text(1.1,.01,'S.O.','fontsize',16','interpreter','latex')
text(1.1,.005,'N.P.','fontsize',16','interpreter','latex')
text(1.1,.0025,'T.P.','fontsize',16','interpreter','latex')
text(1.1,.00125,'N.A.','fontsize',16','interpreter','latex')
text(1.1,.000625,'T.A.','fontsize',16','interpreter','latex')
x = linspace(-.5,2);
f = -0.4543.*x.^2 + 0.2513;
hold on;
lgnd = plot(10.^x,10.^f,'k','linewidth',2)
lgd = legend(lgnd,'$\alpha_{Si} = 1.8e^{-0.2 \ln(Si)^2}$');
set(lgd,'interpreter','latex','fontsize',16)

%{
Parameter fit using cftool, x = log(Si80), y = log(k.^g.*beta) -- 
General model:
     f(x) = a*x^2+b
Coefficients (with 95% confidence bounds):
       a =     -0.2024  (-0.2074, -0.1975)
       b =      0.5971  (0.5581, 0.6361)

Goodness of fit:
  SSE: 7242
  R-square: 0.4941
  Adjusted R-square: 0.494
  RMSE: 1.053
%}