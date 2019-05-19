close all; clear all; clc

%% Constants
%Image operations
alphaUnsharp = 10;
stepSize = 25; %amount of starting point of snake

%snake model
alpha = 1.0;
beta = 0.2;
gamma = -100;
epsilon = 0.0002;
deltaT = 0.2;

%EdgeDetection
useSobel = true; %for CannyFilter set to false
thresHoldVal = 220; %for sobel Filter
    
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

input_sharpened = imageOperators.performSharpening(input,alphaUnsharp);

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

[xVals_opt,yVals_opt,initRadiusSnake,xCenter,yCenter] = snakeHelper.calcInitialSnakeVals(xpol,ypol,stepSize);

% plot optimized center point
hold on,plot(xCenter, yCenter,'g*')

% This are the initial values for the snake, estimated from the user input
hold on, plot(xVals_opt,yVals_opt,'g-')

%% Initialise snake model
% !! xVals are the columns and yVals are the rows in the image!!!
snake = snakeModel.create(alpha,beta,gamma,xVals_opt,yVals_opt, input_medianFil,useSobel,thresHoldVal,initRadiusSnake,xCenter,yCenter);
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
imshow(input)
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






