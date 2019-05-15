close all; clear all; clc

%% Constants
%Image operations
alphaUnsharp = 10;
stepSize = 25;

%snake model
alpha = 0.5;
beta = 0.1;
gamma = 0.2;

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

[x,y] = getline(fig2);
[M,xpol,ypol] = roipoly(input_med,x,y);

%hold on, plot(xpol,ypol)

[xVals_opt,yVals_opt,initRadiusSnake,xCenter,yCenter] = snakeHelper.calcInitialSnakeVals(xpol,ypol,stepSize);

% plot optimized center point
hold on,plot(xCenter, yCenter,'g*')

% This are the initial values for the snake, estimated from the user input
hold on, plot(xVals_opt,yVals_opt,'g-')
%% Initialize snake model
snake = snakeModel.create(alpha,beta,gamma,xVals_opt,yVals_opt, input_med);
image = snake.edgeImage;
figure(4), imshow(snake.edgeImage)


figure(5)
plot3(xVals_opt, yVals_opt, snake.energyValsInit)


initEnergy = snake.totalEnergyInit;
%[snake,gradMag] = snake.minimizeEnergy(0.5);
iterationsteps = 200;
snakeEnergies = zeros(1,iterationsteps);
for i=1:iterationsteps
    
    snake = snake.minimizeEnergy(0.5);
    snakeEnergies(i) = snake.totalEnergy;
    figure(fig2)
    
    if mod(i,2) ~= 0
        color = 'r.';
    else
        color = 'b.';
    end
    xValFinal = snake.finalXVals;
    yValFinal = snake.finalYVals;
    xVal = snake.xVals;
    yVal = snake.yVals;
   
    plot(xValFinal,yValFinal, color)
 
end
endEnergy = snake.totalEnergy;
figure(7)
plot(1:iterationsteps,snakeEnergies)

%GradMagnitude of edge image
%figure(6)
%imshow(gradMag, [])

%figure(fig2)
%plot(snake.yVals,snake.xVals, 'r.')
