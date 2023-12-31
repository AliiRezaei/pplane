clc
clear
close all

%% Define System

sys = @(t,x) [x(2); -2*x(1) - x(2)^3];
dt = 0.01;                       % time step
SimTime = 10;                    % simulation time
t = (0:dt:SimTime)';             % time span
x0 = randn(2,1);                 % initial condition
xLim = [-5 5];                   % x axis limit
yLim = [-5 5];                   % y axis limit
pp = pplane(sys,t,x0,xLim,yLim); % pplane object

%% Phase-Plane Analysis and Plots

figure
pp.plotPhaseTraj();
hold on
pp.isocline(10)
% pp.plotQuiver(10)
xlabel('$x_{1}$','Interpreter','latex','FontSize',14);
ylabel('$x_{2}$','Interpreter','latex','FontSize',14);
title('$Phase Plane$','Interpreter','latex','FontSize',14);
legend('Phase Trajectory','Start Point','End Point');