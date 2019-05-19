close all; clear all; clc

%% Constants
%Image operations
alphaUnsharp = 10;
stepSize = 50; %amount of starting point of snake

%snake model
alpha = 1.0;
beta = 0.2;
gamma = -100;
epsilon = 0.02;
deltaT = 1.0;

%EdgeDetection
useSobel = true; %for CannyFilter set to false
thresHoldVal = 150; %for sobel Filter

%% Input
%input_img = imread('data/simpleObjects.jpg');
%input_img = imread('tire.tif');
%input_img = imread('eight.tif');
%input_img = imread('AT3_1m4_07.tif'); %use the sobel filter; threshold:100
%input_img = imread('hands1.jpg');
input_img = imread('toyobjects.png');
%input_img = imread('pillsetc.png'); %sobel filter; threshold:150 --> good example
subplot(3,3,1)
imshow(input_img)

%% Image operations
[~,~,k]=size(input_img);
if k > 2
    input_grey = imageOperators.convertToGrey(input_img);
else
    input_grey = uint8(input_img);
end

subplot(3,3,2)
imshow(input_grey)

input_med = imageOperators.medianFilter(input_grey);


%% Initialse snake
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
%% Initialise snake model
% !! xVals are the columns and yVals are the rows in the image!!!
snake = snakeModel.create(alpha,beta,gamma,xVals_opt,yVals_opt, input_med,useSobel,thresHoldVal);
image = snake.edgeImage;

figure(4), imshow(snake.edgeImage)

snakeEnergImg = snake.energyImage;


%% Minimize snake energy
initEnergy = snake.totalEnergyInit;
iterationsteps = 1;
snakeEnergies = [];
snakeXYVals = {};
figure(fig2)
color = 'b-';
snakeCont = plot(snake.xVals,snake.yVals, color);
prevEnergy = snake.totalEnergy;
runLoop = true;
cellIndex = 1;
while runLoop 
   
    snake = snake.minimizeEnergy(deltaT);
    snakeEnergies(iterationsteps) = snake.totalEnergy;
    if mod(iterationsteps,10) == 0
        snakeXYVals{1, cellIndex} = snake.xVals;
        snakeXYVals{2, cellIndex} = snake.yVals;
        cellIndex = cellIndex + 1;
    end
    
    iterationsteps = iterationsteps + 1;
    actualEnergy = snake.totalEnergy;
    if abs(prevEnergy-actualEnergy) > epsilon
        runLoop = true;
        prevEnergy = actualEnergy;  
    else
        runLoop = false;
    end
    
    if mod(iterationsteps,2) ~= 0
        color = 'r-';
    else
        color = 'b-';
    end
    figure(fig2)
    delete(snakeCont)
    snakeCont = plot(snake.xVals,snake.yVals, color);
    
    
end

endEnergy = snake.totalEnergy;
figure(7)
plot(1:length(snakeEnergies),snakeEnergies)

fig8 = figure(8);
imshow(input_grey)
figure(fig8)
hold on,plot(snake.xVals,snake.yVals,'g-')
title(["Final snake after ",num2str(iterationsteps), "iterations"])

figure(9)
[rows,columns] = size(snakeXYVals);
n = length(snakeXYVals{1,1});
 for i=1:columns-1
     X = snakeXYVals{1,i};
     Y = snakeXYVals{2,i};
     X_next = snakeXYVals{1,i+1};
     Y_next = snakeXYVals{2,i+1};
     diffX = X_next - X;
     diffY = Y_next - Y;
     plot(X,Y),hold on
     quiver(X,Y,diffX,diffY,0,'color',[1 0 0]),hold on
end



