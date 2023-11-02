classdef pplane < handle
    % The Phase Plane Analysis and Plots Class is a MATLAB class ...
    % designed to facilitate the analysis of second order linear ...
    % or non-linear dynamical systems. By passing the system equations ...
    % as a handle function to the pplane and defining the time span, ...
    % initial condition, and axis limitations, users can utilize the ...
    % class methods to plot phase trajectories and vector fields using ...
    % the isocline method. This class provides a user-friendly interface ...
    % for exploring the behavior of second order systems and can be a ..
    % valuable tool for non-linear systems researchers.
    % author : Ali Rezaei
    % date : November, 1, 2023
    % MATLAB version : R2023a

    properties
        sys;
        tspan;
        InitCond;
        sol = struct('States',[],'Time',[]);
        x_lim;
        y_lim;
    end % end of public properties

    properties(GetAccess = private , SetAccess = private)
        symbolic_sys;
        symbolic_states;
    end % end fo private properties

    properties(Constant,Hidden)
        sys_order = 2;
    end % end of constant properties

    methods
        function obj = pplane(sys,t,x0,x_lim,y_lim)
            % Object Constructor.
            % Inputs : 
            %   sys   --> a second order system (must be a handle function)
            %   tspan --> time span
            %   x0    --> initial condition
            % Outputs : pplane object

            % check requirements :
            if ~isa(sys,'function_handle')
                msg = sprintf(['Error using pplane \n ' ...
                    'please pass the system equations as a handle function']);
                error(msg);
            end
            if numel(x0) ~= 2
                msg = sprintf(['Error using pplane \n ' ...
                   'system or initial condition dimension not correct']);
                error(msg);
            end

            obj.sys = sys;
            obj.symbolic_states = sym('x',[obj.sys_order , 1],'real');
            obj.symbolic_sys = obj.sys([],obj.symbolic_states);
            obj.tspan = t;
            obj.InitCond = x0;
            obj.x_lim = x_lim;
            obj.y_lim = y_lim;
            
            [T,X] = obj.solveSys();
            obj.sol.States = X;
            obj.sol.Time = T;
        end % end of constructor

        function [T,X] = solveSys(obj)
            % This method allows you to solve your system equation.
            % Inputs  : pplane object
            % Outputs :
            %    X --> numercal solution of system
            %    T --> time interval

            nt = numel(obj.tspan); % number of element of t
            dt = round(max(obj.tspan) - min(obj.tspan)) / nt; % time resolution
            [T,X] = obj.RK4(obj.tspan,dt,obj.InitCond); % solve system equation using RK4
        end % end of solveSys method

        function isocline(obj,n,varargin)
            % This method is the implementation of the isocline method ...
            % for plotting vector fields of any second-order system.
            % Inputs :
            %   n --> number of arrows in each axis
            % Outputs : Nothing

            % check isoclines show status :
            if ismember(varargin,'ShowIsoclines')
                ShowIsoclines = true;
            else
                ShowIsoclines = false;
            end

            % initialize x , y axis :
            xAxis = linspace(obj.x_lim(1),obj.x_lim(2),n);
            yAxis = linspace(obj.y_lim(1),obj.y_lim(2),n);
            
            % main loop :
            for i = xAxis
                for j = yAxis
                    f = double(subs(obj.sys([],[i; j])));
                    % determine slope :
                    if f(1) == 0
                        slope = sign(f(2)) * 90;
                    elseif f(2) == 0 && sign(f(1)) > 0
                        slope = 0;
                    elseif f(2) == 0 && sign(f(1)) < 0
                        slope = 180;
                    else
                        slope = atan2(f(2),f(1)) * 180 / pi;
                    end
                    obj.plotArrow(i,j,slope,'Color','r', ...
                        'HorizontalAlignment','Center');
                    if ShowIsoclines
                        slope = tan(slope * pi / 180); % degree to radian
                        % solve f(2) = alpha * f(1) :
                        curve = vpasolve(obj.symbolic_sys(2) == slope * obj.symbolic_sys(1),obj.symbolic_states(2),4);
                        hold on
                        h = fplot(curve,obj.x_lim);
                        h.Color = zeros(1,3);
                        h.LineWidth = 0.2;
                        hold off
                    end
                end % end of internal loop
            end % end of main loop

            % set limitation of current axis :
            currentAxis = gca;
            currentAxis.XLim = obj.x_lim;
            currentAxis.YLim = obj.y_lim;
        end % end of isocline method

        function plotQuiver(obj,n)
            % This method plot vector fields.
            % using MATLAB quiver command.
            % Inputs :
            %   n  --> number of arrows in each axis
            % Outputs : Nothing

            % initialize vector fields :
            x1 = linspace(obj.x_lim(1),obj.x_lim(2),n);
            x2 = linspace(obj.y_lim(1),obj.y_lim(2),n);
            [x1,x2] = meshgrid(x1,x2);
            u = zeros(size(x1));
            v = zeros(size(x2));

            % calculate velocities :
            for i = 1:numel(x1)
                xprim = obj.sys(0,[x1(i);x2(i)]);
                u(i) = xprim(1);
                v(i) = xprim(2);
            end
            hold on
            quiver(gca,x1,x2,u,v,'r') % plot vector fields
            hold off
        end % end of plotQuiver method

        function plotPhaseTraj(obj)
            % Phase Trajectory Plotter.
            hold on
            plot(gca,obj.sol.States(:,1),obj.sol.States(:,2),...
                'Color',[0 0.447 0.741],'LineWidth',1.5)
            StartPoint = obj.InitCond;
            EndPoint = [obj.sol.States(end,1),obj.sol.States(end,2)];
            plot(gca,StartPoint(1),StartPoint(2),'ro',...
                'MarkerFaceColor','r')
            plot(gca,EndPoint(1),EndPoint(2),'go',...
                'MarkerFaceColor','g')
            hold off
        end % end of plotPhaseTraj method
    end % end of public methods

    methods(Access = private)
        function [t,Y] = RK4(obj,tSpan,dt,Y0)
            % This method is the Implementation of the 4th-Order Runge-Kutta Method
            % for the Numeric Solve of Differetial Equations...
            % Inputs:
            %     tSpan --> time span for the ODE solve
            %     dt    --> time resolution
            %     x0    --> initial condition
            %
            % Outputs:
            %     T --> simulation time
            %     Y --> ODE numerical solution
            % copyright :
            %     Bnaan Kiamanesh, (https://github.com/BanaanKiamanesh)

            StepNum  = floor((max(tSpan) - min(tSpan)) / dt) + 1;    % number of simulation steps
            StateNum = numel(Y0);                                    % number of sys states
            odeFun = obj.sys; % ode function

            % initialize solution and time vars
            Y = zeros(StateNum, StepNum);
            t = (min(tSpan) : dt : max(tSpan))';
            Y(:, 1) = Y0;

            % main loop
            for k = 1:StepNum - 1
                K1 = dt * feval(odeFun, t(k), Y(:, k));
                K2 = dt * feval(odeFun, t(k) + dt/2 , Y(:, k) + K1 / 2);
                K3 = dt * feval(odeFun, t(k) + dt/2 , Y(:, k) + K2 / 2);
                K4 = dt * feval(odeFun, t(k) + dt   , Y(:, k) + K3);
                Y(:, k+1) = Y(:, k) + (K1 + K2 * 2 + K3 * 2 + K4) / 6;
            end
            % transpose to make this and ODE45 look alike
            Y = transpose(Y);
        end % end of RK4 method

        function plotArrow(~,x,y,deg,varargin)
            % This method plot vector feild arrows
            % Inputs :
            %   x   --> x position of arrow
            %   y   --> y position of arrow
            %   deg --> orientation of arrow (degree)

            h = text(x,y,'\rightarrow',varargin{:});
            set(h,'Rotation',deg);
        end % end of plotArrow method
    end % end of private methods
end % end of class