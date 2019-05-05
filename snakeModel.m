classdef snakeModel

    
    methods(Static)
        function [potVal,image] = imageForces(imageData)
                %IMAGEFORCES draws Potential towards an edge of image
                 %   edge detection: x -> x + ddI(x) + dI(x)

            %% Blurring
            %blurred_Image = imageOperators.gaussianFilter(imageData);
            blurred_Image = imageOperators.anisotropicFilter(imageData);
            figure(3)
            subplot(2,3,1)
            imshow(blurred_Image)
            title("Blurred Image")

            %% Edge Detection
            %edge_Image = imageOperators.detectEdgesLaplace(blurred_Image);
            edge_Image = imageOperators.detectEdgesSobel(blurred_Image);
            subplot(2,3,2)
            imshow(uint8(edge_Image))
            title('Edges after Sobel (2nd order)')

            %edge_Image = imageOperators.detectEdgesXY(edge_Image);
            subplot(2,3,3)
            imshow(edge_Image)
            title('final edge detection (1.2nd Order xy, 2.1st Order xy')

            canny_Image = imageOperators.cannyFilter(blurred_Image);
            subplot(2,3,5)
            imshow(canny_Image)
            title('canny filtered')

            %% Inverting
            image = imcomplement(edge_Image); %if snake hits an edge potential gets zero
            subplot(2,3,4)
            imshow(image)
            title('inverted Image')

            potVal = norm(double(image)); %I think this must be the potvalue of the snake and not the image

        end
        
        function [xData,yData] = internalForces(xData,yData,alpha, beta)
            %INTERNALFORCES calculation of tension and stiffness
            %   calculate first and second derivative of snake function (C)
            %   E_intern = (alpha/2) * C' + (beta/2) * C''

            %ToDo: Implement

        end
        
        function [xData,yData] = snakeSegmentation(alpha, beta,xData,yData,imageData)
            %SNAKESEGMENTATION total Snake Energy function
            %   Calculate Coordinates for snake with lowest energy potentials
            %   E(snake) = E_intern + E_image

            %TODO: implement

        end

    end
end

