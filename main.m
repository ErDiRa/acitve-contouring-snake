close all; clear all; clc

%% Constants
%Image operations
alphaUnsharp = 10;

%snake model
alpha = 1.0;
beta = 0.00002;
gamma = -10;

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

hold on, plot(xpol,ypol)

% estimate centerpoint and radius of given user input
[xCenter,yCenter] = snakeHelper.calcCenterOfPoints(xpol,ypol);
plot(xCenter,yCenter,'r.'), hold on

radius = snakeHelper.estimateRadius(xCenter,yCenter,xpol,ypol);

% Optimise circle function given the polygon input of user
[optCircleVals, R] = snakeHelper.optimiseCircleParams(xCenter,yCenter,radius,xpol,ypol);

% plot optimized center point
plot(optCircleVals(1), optCircleVals(2),'g*'), hold on

% plot optimzed circel of user input
stepSize = 50;
% This are the initial values for the snake, estimated from the user input
[xVals_opt, yVals_opt] = snakeHelper.calcCirclePlotVals(optCircleVals(1), optCircleVals(2), optCircleVals(3),stepSize);
plot(xVals_opt,yVals_opt,'g-')

%% Initialize snake model
% !! xVals are the columns and yVals are the rows in the image!!!
snake = snakeModel.create(alpha,beta,gamma,xVals_opt,yVals_opt, input_medianFil,useSobel,thresHoldVal);
image = snake.edgeImage;
[rows,columns] = size(input_medianFil);
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







