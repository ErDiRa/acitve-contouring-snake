classdef snakeModel
    
    properties
        alpha
        beta
        gamma
        xVals
        yVals
        energyVal
    end
    
    methods(Static)
        
        function snake = create(alpha,beta,gamma,xVals,yVals)
            snake = snakeModel;
            snake.alpha = alpha;
            snake.beta = beta;
            snake.gamma = gamma;
            snake.xVals = xVals;
            snake.yVals = yVals;
         end
        
        function edgeImage = prepareEdgeImage(imageData)
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
            edgeImage = imcomplement(edge_Image); %if snake hits an edge potential gets zero
            subplot(2,3,4)
            imshow(edgeImage)
            title('inverted Image')
        end
        
        function imageEnergy = imageForces(dataset,edgeImage)
          
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
    
    methods
        
         
        
         function tension = calcTension(this)
            %ToDO: i need a total tension and a tension for each point
            x_diff = diff(this.xVals)';
            y_diff = diff(this.yVals)';
            norm_XY_Diff = sqrt(x_diff.^2 + y_diff.^2); %ToDO: something wrong
            
            tension = this.alpha * sum(norm_XY_Diff);
         
         end
         
         function stiffness = calcStiffness(this)
            stiffness = 0;
         end
        
    end
end

