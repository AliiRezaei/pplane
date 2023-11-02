clc
clear
close all

%% Define System

epsilon = .1;
sys = @(t,x) [x(2); -x(1)+epsilon*(x(2)-x(2)^3/3)]; % rayliegh's equation
dt = 0.01;                       % time step
SimTime = 50;                    % simulation time
t = (0:dt:SimTime)';             % time span
x0 = randn(2,1);                 % initial condition
xLim = [-3 3];                   % x axis limit
yLim = [-3 3];                   % y axis limit
pp = pplane(sys,t,x0,xLim,yLim); % pplane object

%% Phase-Plane Analysis and PLots

figure
pp.plotPhaseTraj();
hold on
pp.isocline(10);
% pp.plotQuiver(10);
xlabel('$x_{1}$','Interpreter','latex','FontSize',14);
ylabel('$x_{2}$','Interpreter','latex','FontSize',14);
title('$Phase Plane$','Interpreter','latex','FontSize',14);
legend('Phase Trajectory','Start Point','End Point');
