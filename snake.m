close all; clear all; clc
%% Read input image
I = dicomread('./data/ElasticRadExampleData/BrainX/20061201/IM-0001-0009.dcm');

input = imageOperators.convertGreyValsToInt8(I);

figure(1)
subplot(2,2,1)
imshow(input)
title('Original Image')
%imtool(input)

%figure
%imhist(input,65)
%title('Histogram')
% 

%Unsharpmasking
alpha = 10;
input_sharpened = imageOperators.performSharpening(input,alpha);

subplot(2,2,2)
imshow(input_sharpened)
title('Snake Input image sharpened')

%Windowing
gmin = 50;
gmax = 120;
maxPixelVal = 255;
%input_eq = adapthisteq(input_sharpened); %Histogram equalisation CLAHE
input_win = imageOperators.performWindowing(input_sharpened,gmin,gmax,maxPixelVal);
% figure
% imhist(input_win,65)
% title('Histogram windowed')
subplot(2,2,3)
imshow(input_win)
title_str = ['Snake Input image [',num2str(gmin),',',num2str(gmax),']'];
title(title_str)

%Medianfilter (reduce noise) (will preserve edges!!)
input_medianFil = imageOperators.medianFilter(input_win);

subplot(2,2,4)
imshow(input_medianFil)
title('Snake Input Image Median Filter')

fig2 = figure(2);
imshow(input_medianFil)
title('Snake Input Image Median Filter')

%% Create initial snake
[x,y] = getline(fig2);
[M,xpol,ypol] = roipoly(input_medianFil,x,y);

hold on, plot(xpol,ypol)

% estimate centerpoint and radius of given user input
[xCenter,yCenter] = calcCenterOfPoints(xpol,ypol);
plot(xCenter,yCenter,'r.'), hold on

radius = estimateRadius(xCenter,yCenter,xpol,ypol);

% Optimise circle function given the polygon input of user
[optCircleVals, R] = optimiseCircleParams(xCenter,yCenter,radius,xpol,ypol);

% plot optimized center point
plot(optCircleVals(1), optCircleVals(2),'g*'), hold on

% plot optimzed circel of user input
stepSize = 50;
% This are the initial values for the snake, estimated from the user input
[xVals_opt, yVals_opt] = calcCirclePlotVals(optCircleVals(1), optCircleVals(2), optCircleVals(3),stepSize);
plot(xVals_opt,yVals_opt,'g-')

%TODO: use the snake model to minimize the energy of the snake while
%changing the x and y values

%% Smooth image and detect edges (inverted)
% Smooth image and detect edges
[potVal, image_edge] = snakeModel.imageForces(input_medianFil);
figure(4), imshow(image_edge)

%% Functions
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
    circle_fun = @(inputs) circleFun(inputs(1),inputs(2), inputs(3), xData,yData);
    opt_params = [x0,y0,r];
    [optVals, R] = lsqnonlin(circle_fun, opt_params);
end

% TODO: this is also the function which has to be derived C = (x(s), y(s))
% what does not actually matter be cause we use numerical differentiation
function [xVals, yVals] = calcCirclePlotVals(x0,y0, r, steps)
    %using parameter form
    s = 0:pi/steps:2*pi;
    
    xVals = r * cos(s) + x0;
    yVals = r * sin(s) + y0;
    
end

