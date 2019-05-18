close all; clear all; clc

%% Constants
%Image operations
alphaUnsharp = 10;
stepSize = 25;

%snake model
alpha = 1.0;
beta = 0.2;
gamma = -100;

%EdgeDetection
useSobel = true; %for CannyFilter set to false
thresHoldVal = 150; %for sobel Filter

%% Input
input_img = imread('data/simpleObjects.jpg');
%input_img = imread('tire.tif');
%input_img = imread('eight.tif');
%input_img = imread('AT3_1m4_07.tif'); %use the sobel filter; threshold:100
%input_img = imread('hands1.jpg');
%input_img = imread('toyobjects.png');
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
% !! xVals are the columns and yVals are the rows in the image!!!
snake = snakeModel.create(alpha,beta,gamma,xVals_opt,yVals_opt, input_med,useSobel,thresHoldVal);
image = snake.edgeImage;
[rows,columns] = size(input_med);
figure(4), imshow(snake.edgeImage)

snakeEnergImg = snake.energyImage;

figure(5)
mesh(1:columns,1:rows,snakeEnergImg)

figure(6)
plot3(xVals_opt, yVals_opt, snake.energyValsInit)


initEnergy = snake.totalEnergyInit;
iterationsteps = 500;
snakeEnergies = zeros(1,iterationsteps);
figure(fig2)
color = 'b-';
snakeCont = plot(snake.xVals,snake.yVals, color);
for i=1:iterationsteps
    
    snake = snake.minimizeEnergy(0.2);
    snakeEnergies(i) = snake.totalEnergy;
    
    
    if mod(i,2) ~= 0
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
plot(1:iterationsteps,snakeEnergies)
