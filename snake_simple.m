close all; clear all; clc

%% Constants
%Image operations
alphaUnsharp = 10;
stepSize = 50;

%snake model
alpha = 1.0;
beta = 1.0;
gamma = -1.0;

%% Input
input_img = imread('data/simpleObjects.jpg');

subplot(3,3,1)
imshow(input_img)

%% Image operations
input_grey = imageOperators.convertToGrey(input_img);
subplot(3,3,2)
imshow(input_grey)

input_med = imageOperators.medianFilter(input_grey);


%% Init snake
fig2 = figure(2);

imshow(input_med)
[height,width]= size(input_med);

[x,y] = getline(fig2);
[M,xpol,ypol] = roipoly(input_med,x,y);

%hold on, plot(xpol,ypol)

[xVals_opt,yVals_opt,xCenter,yCenter] = snakeHelper.calcInitialSnakeVals(xpol,ypol,stepSize);

% plot optimized center point
hold on,plot(xCenter, yCenter,'g*')

% This are the initial values for the snake, estimated from the user input
hold on, plot(xVals_opt,yVals_opt,'g-')
%% Initialize snake model
%xvals and yvals has to be swapped: x==height, y==width 
%--> in a plot x expresses the width and y the height
snake = snakeModel.create(alpha,beta,gamma,yVals_opt,xVals_opt, input_med);

figure(4), imshow(snake.edgeImage)
[tensionVals, stiffnessVals,potVals, energyVals, totalEnergyInit] = snake.calcEnergyVals();
[secDerivX, secDerivY, fourthDerivX, fourthDerivY,imageGradVals,gradMag] = snake.minimizeEnergy(1.0);

figure(5)
plot3(xVals_opt, yVals_opt, energyVals)

figure(6)
imshow(gradMag, [])
