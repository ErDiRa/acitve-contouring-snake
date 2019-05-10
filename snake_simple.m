close all; clear all; clc
input_img = imread('data/simpleObjects.jpg');

subplot(3,3,1)
imshow(input_img)

input_grey = imageOperators.convertToGrey(input_img);
subplot(3,3,2)
imshow(input_grey)

thres = 200;
input_thres = imageOperators.performThresholding(input_grey,thres);
subplot(3,3,3)
imshow(input_thres)

[potVal, image_edge] = snakeModel.imageForces(input_thres);
figure, imshow(image_edge)

image_edgeBlurred = imageOperators.gaussianFilter(image_edge);
figure, imshow(image_edgeBlurred)

