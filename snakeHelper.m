classdef snakeHelper
    %SNAKEHELPER contains functions to help initializing snake
   
    methods(Static)
        
        function [xVals, yVals,xCenter,yCenter] = calcInitialSnakeVals(xpol,ypol,stepSize)
            % estimate centerpoint and radius of given user input
            [xCenterPol,yCenterPol] = snakeHelper.calcCenterOfPoints(xpol,ypol);

            radius = snakeHelper.estimateRadius(xCenterPol,yCenterPol,xpol,ypol);
            
            % Optimise circle function given the polygon input of user
            optCircleVals = snakeHelper.optimiseCircleParams(xCenterPol,yCenterPol,radius,xpol,ypol);

            % optimized center point
            xCenter = optCircleVals(1);
            yCenter = optCircleVals(2);
      
            % This are the initial values for the snake, estimated from the user input
            [xVals, yVals] = snakeHelper.calcCirclePlotVals(optCircleVals(1), optCircleVals(2), optCircleVals(3),stepSize);
         
        end
      
        function [xCenter,yCenter] = calcCenterOfPoints(xData,yData)
            %estimates center for fit of circular function
            xCenter = mean(xData);
            yCenter = mean(yData);
        end

        function radius = estimateRadius(xCenter,yCenter,xData,yData)
            % estimates radius for fit of circular function
            distances = zeros(length(xData),1);

            for i=1:length(xData)
                d = sqrt( (xData(i) - xCenter)^2 + (yData(i)-yCenter)^2);
                distances(i) = d;
            end

            radius = mean(d);

        end

        function f_circle = circleFun(a, b, r, xData, yData)
            % https://en.wikipedia.org/wiki/Circle#Equations
            n = length(xData);
            f_circle = zeros(n,1);
            for i=1:n
                f_circle(i) = abs((xData(i)-a)^2 + (yData(i)-b)^2 - r^2);
            end
        end

        function [optVals, R] = optimiseCircleParams(x0, y0, r, xData, yData)
            %Define inputs for least square function of matlab
            %Returns opt values for: center points of circle and radius
            circle_fun = @(inputs) snakeHelper.circleFun(inputs(1),inputs(2), inputs(3), xData,yData);
            opt_params = [x0,y0,r];
            [optVals, R] = lsqnonlin(circle_fun, opt_params);
        end

        function [xVals, yVals] = calcCirclePlotVals(x0,y0, r, steps)
            %using parameter form
            s = 0:pi/steps:2*pi;

            xVals = round(r * cos(s) + x0);
            yVals = round(r * sin(s) + y0);

        end
    end
end

