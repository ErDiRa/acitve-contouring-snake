classdef snakeModel
    % Author: Stefan Karlin
    % link for numerical differentiation: 
    % https://www.sphenisc.com/doku.php/software/algo/finitedifference
    % used main literature:
    % EVERYTHING YOU ALWAYS WANTED TO KNOW ABOUT SNAKES 
    % (BUT WERE AFRAID TO ASK) Jim Ivins  & John Porrill (in pdf dir)
    properties
        %Snake parameters
        alpha
        beta
        gamma
        %Contour coordinates
        xVals
        yVals
        %eEnergy vals of snake
        energyVals
        totalEnergy
        energyValsInit
        totalEnergyInit
        edgeImage
        energyImage
        %Values of init snake to restrict further calc points
        radius
        xCenter
        yCenter
    end
    
    methods(Static)
        
        function snake = create(alpha,beta,gamma,xVals,yVals,imageData,useSobel,thresVal,radius,xCenter,yCenter)
            % Initialises snake model
            snake = snakeModel;
            snake.alpha = alpha;
            snake.beta = beta;
            snake.gamma = gamma;
            snake.xVals = xVals;
            snake.yVals = yVals;
            snake.radius = radius;
            snake.xCenter = xCenter;
            snake.yCenter = yCenter;
            snake.edgeImage = snakeModel.prepareEdgeImage(imageData,useSobel,thresVal);
            [snake.energyValsInit, snake.totalEnergyInit ] = ...
                snakeModel.calcEnergyVals(xVals,yVals,snake.edgeImage, snake.alpha,snake.beta,snake.gamma);
            snake.totalEnergy = snake.totalEnergyInit;
            snake.energyVals = snake.energyValsInit;
        end
        
    end
    
    methods(Access = private, Static)
        % Block contains functions to calc all needed values for the snake
        
        function edgeImage = prepareEdgeImage(imageData,useSobel,thresVal)
            % function uses either Sobel or Canny Filter to detect the
            % edges
            % Canny Filter needs more optimisation to ignore unnecessary
            % edges (thresholding)
            %% Blurring
            %blurred_Image = imageOperators.gaussianFilter(imageData);
            % use a anisotropicFilter to conserve edges (better then gauss)
            blurred_Image = imageOperators.anisotropicFilter(imageData);
            figure(3)
            subplot(2,3,1)
            imshow(blurred_Image,[])
            title("Blurred Image")

            %% Sobel Edge Detection
            edge_Image = imageOperators.detectEdgesSobel(blurred_Image);
            subplot(2,3,2)
            imshow(uint8(edge_Image))
            title('Edges after Sobel (2nd order)')
            
            edge_Image = imageOperators.performThresholding(edge_Image,thresVal);
            subplot(2,3,3)
            imshow(uint8(edge_Image))
            title('Edges after Sobel and Thres')
            
            %% Canny Edge Detection
            canny_Image = imageOperators.cannyFilter(blurred_Image);
            canny_Image_blurred = imageOperators.anisotropicFilter(canny_Image);
            subplot(2,3,4)
            imshow(canny_Image)
            title('canny filtered and blurred')

            %% Inverting
            if useSobel
                edgeImage = edge_Image;
            else
                edgeImage = canny_Image_blurred;
            end
           
            subplot(2,3,5)
            imshow(edgeImage)
            title('final edge Image')
        
         end 
        
        function [gradMagValues,gradMagDirection,gradMag] = calcGradientMagnitude(edgeImage,xVals,yVals)
            % calcs gradient magnitude of edge image
            % needed to calc gradient vals
            [gradMag, gradMagDir] = imgradient(edgeImage);
            gradMagValues = zeros(1,length(xVals));
            gradMagDirection = zeros(1, length(yVals));
            for i=1:length(xVals)
                xpos = round(xVals(i));
                ypos = round(yVals(i));
                %!! xVals are columns and yVals are the rows !!, due to
                %getline() function use in snake_simple.m or snake_brain.m
                gradMagValues(i) = gradMag(ypos,xpos);
                gradMagDirection = gradMagDir(ypos,xpos);
            end
            
            gradMag = uint8(gradMag);
            gradMagValues = uint8(gradMagValues);
            
        end
        
        function [secDerivX, secDerivY] = calcSecondDerivative(xVals,yVals)
             % calcs second derivatives of  x and y vals of snake
             % needed to calc gradient vals
             n = length(xVals);
             secDerivX = zeros(1,n);
             secDerivY = zeros(1,n);
             for i = 1:n
                  % needs Index checking because we have a circular closes
                 % snake
                 if i == 1
                    v_prev = [xVals(n-1); yVals(n-1)]; 
                    v = [xVals(i); yVals(i)];
                    v_next = [xVals(i+1);yVals(i+1)];

                    %2nd derivative
                    secDerivX(i) = v_prev(1) - 2*v(1) + v_next(1);
                    secDerivY(i) = v_prev(2) - 2*v(2) + v_next(2);
                 
                 elseif i < n
                    
                    v_prev = [xVals(i-1); yVals(i-1)];
                    v = [xVals(i); yVals(i)];
                    v_next = [xVals(i+1);yVals(i+1)];

                    %2nd derivative
                    secDerivX(i) = v_prev(1) - 2*v(1) + v_next(1); 
                    secDerivY(i) = v_prev(2) - 2*v(2) + v_next(2);
                elseif i == n
                     secDerivX(i) = secDerivX(1);
                     secDerivY(i) = secDerivY(1);
                 end
   
             end
         end
         
        function [fourthDerivX, fourthDerivY] = calcFourthDerivative(xVals,yVals)
             % calcs fourth derivatitves of x and y vals of snake 
             % needed to calc gradient vals
             n = length(xVals);
             fourthDerivX = zeros(1,n);
             fourthDerivY = zeros(1,n);
             for i = 1:n
                 % needs Index checking because we have a circular closes
                 % snake
                 if i == 1
                    v_pprev = [xVals(n-2); yVals(n-2)];
                    v_prev = [xVals(n-1); yVals(n-1)]; 
                    v = [xVals(i); yVals(i)];
                    v_next = [xVals(i+1);yVals(i+1)];
                    v_nnext = [xVals(i+2);yVals(i+2)];

                    %4th derivative
                    fourthDerivX(i) = v_pprev(1) - 4* v_prev(1) + 6 * v(1) - 4* v_next(1) + v_nnext(1);
                    fourthDerivY(i) = v_pprev(2) - 4* v_prev(2) + 6 * v(2) - 4* v_next(2) + v_nnext(2);
                    
                 elseif i == 2
                    v_pprev = [xVals(n-1); yVals(n-1)];
                    v_prev = [xVals(1); yVals(1)]; 
                    v = [xVals(i); yVals(i)];
                    v_next = [xVals(i+1);yVals(i+1)];
                    v_nnext = [xVals(i+2);yVals(i+2)];

                    %4th derivative
                    fourthDerivX(i) = (v_pprev(1) - 4* v_prev(1) + 6 * v(1) - 4* v_next(1) + v_nnext(1));
                    fourthDerivY(i) = v_pprev(2) - 4* v_prev(2) + 6 * v(2) - 4* v_next(2) + v_nnext(2);     
                    
                 elseif i < n - 2
                    v_pprev = [xVals(i-2); yVals(i-2)];
                    v_prev = [xVals(i-1); yVals(i-1)];
                    v = [xVals(i); yVals(i)];
                    v_next = [xVals(i+1);yVals(i+1)];
                    v_nnext = [xVals(i+2);yVals(i+2)];
                    
                    %4th derivative
                    fourthDerivX(i) = v_pprev(1) - 4* v_prev(1) + 6 * v(1) - 4* v_next(1) + v_nnext(1);
                    fourthDerivY(i) = v_pprev(2) - 4* v_prev(2) + 6 * v(2) - 4* v_next(2) + v_nnext(2);
                           
                 elseif i == n -1 
                     v_pprev = [xVals(i-2); yVals(i-2)];
                     v_prev = [xVals(i-1); yVals(i-1)];
                     v = [xVals(i); yVals(i)];
                     v_next = [xVals(1);yVals(1)];
                     v_nnext = [xVals(2);yVals(2)];
                     
                     %4th derivative
                     fourthDerivX(i) = v_pprev(1) - 4* v_prev(1) + 6 * v(1) - 4* v_next(1) + v_nnext(1);
                     fourthDerivY(i) = v_pprev(2) - 4* v_prev(2) + 6 * v(2) - 4* v_next(2) + v_nnext(2);
                
                 elseif i == n
                    fourthDerivX(i) = fourthDerivX(1);
                    fourthDerivY(i) = fourthDerivY(1);               
                 end
                 
                
   
             end
        end
         
        function [tensionValues,totalTension] = calcTension(xVals,yVals)
            
            n = length(xVals);
            tensionValues = zeros(1,n);
                
            for i=1:n
                 % needs index checks: because we have a circular closed
                 % snake
                if i == 1
                    v_prev = [xVals(n-1); yVals(n-1)];
                    v = [xVals(i); yVals(i)];
                    v_next = [xVals(i+1);yVals(i+1)];
                    
                    deriveX = - 0.5*v_prev(1) + 0*v(1) + 0.5*v_next(1);
                    deriveY = - 0.5* v_prev(2) + 0*v(2) + 0.5*v_next(2);

                    tensionValues(i) = (sqrt(deriveX^2 + deriveY^2))^2;
               
                elseif i < n
                    v_prev = [xVals(i-1); yVals(i-1)];
                    v = [xVals(i); yVals(i)];
                    v_next = [xVals(i+1);yVals(i+1)];
                    %1st derivative (forward difference)
                    deriveX = - 0.5*v_prev(1) + 0*v(1) + 0.5*v_next(1);
                    deriveY = - 0.5* v_prev(2) + 0*v(2) + 0.5*v_next(2);

                    tensionValues(i) = (sqrt(deriveX^2 + deriveY^2))^2;
                elseif i == n
                    tensionValues(i) = tensionValues(1);
                end
                
            end
            
            totalTension = sum(tensionValues);
         
         end
         
        function [stiffnessVals,stiffnessTotal] = calcStiffness(xVals,yVals)
             n = length(xVals);
             stiffnessVals = zeros(1, n);
             
             for i = 1:n
                 % needs Index checks: because we have a circular closed
                 % snake
                 if i == 1
                    v_prev = [xVals(n-1); yVals(n-1)]; %because at index n is the same value as 1 (close circle)
                    v = [xVals(i); yVals(i)];
                    v_next = [xVals(i+1);yVals(i+1)];

                    %2nd derivative
                    deriveX = v_prev(1)  - 2*v(1) +  v_next(1);
                    deriveY = v_prev(2) - 2*v(2) + v_next(2);

                    stiffnessVals(i) = (sqrt(deriveX^2 + deriveY^2))^2;
                 
                 elseif i < n 
                    
                    v_prev = [xVals(i-1); yVals(i-1)];
                    v = [xVals(i); yVals(i)];
                    v_next = [xVals(i+1);yVals(i+1)];

                    %2nd derivative
                    deriveX = v_prev(1) - 2*v(1) + v_next(1);
                    deriveY = v_prev(2) - 2*v(2) + v_next(2);

                    stiffnessVals(i) = (sqrt(deriveX^2 + deriveY^2))^2;
                
                 else
                    stiffnessVals(i) = stiffnessVals(1);
                    
                end
             end
            
            stiffnessTotal = sum(stiffnessVals);
         end
         
        function [imageEnergyVals,imageEnergyTotal] = calcImageForces(xVals,yVals,edgeImage)
            % return image values of edge detected image
            imageEnergyVals = zeros(1, length(xVals));
            for i=1:length(xVals)
                xpos = round(xVals(i));
                ypos = round(yVals(i));
                %!! xVals are columns and yVals are the rows !!, due to
                %getline() function use in snake_simple.m or snake_brain.m
                imageEnergyVals(i) = double((edgeImage(ypos,xpos)))^2;
                 
            end
            
            imageEnergyTotal = sum(imageEnergyVals);
        end
        
        function [gradientEnergyX,gradientEnergyY] = calcGradientVals(secDerivX,secDerivY, fourthDerivX,...
                fourthDerivY,imageGradVals,alpha,beta,gamma)
            % calculates the gradient values for x and y
            n = length(secDerivX);
            gradientEnergyX = zeros(1,n);
            gradientEnergyY = zeros(1,n);
            
            for i=1:n
                gradientEnergyX(i) = alpha*secDerivX(i) - beta*fourthDerivX(i) - gamma*double(imageGradVals(i));
                gradientEnergyY(i) = alpha*secDerivY(i) - beta*fourthDerivY(i) - gamma*double(imageGradVals(i));
            end
            
            
        end
        
        
        function [energyVals,totalEnergy] = calcEnergyVals(xVals,yVals,...
                edgeImage,alpha,beta,gamma)
             % function calcs total energy and energy for each points of
             % the snake
             [tensionValsNew, ~] = snakeModel.calcTension(xVals,yVals);

             n = length(tensionValsNew(1,:));
           
             [stiffnessValsNew, ~] = snakeModel.calcStiffness(xVals, yVals);

             
             [edgePotentialsNew,~] = snakeModel.calcImageForces(xVals,yVals,edgeImage);
  
             
             energyVals = zeros(1,n);
             totalEnergy = 0;
             for i=1:n
                 energyValNew = (alpha/2)*tensionValsNew(i) + (beta/2)*stiffnessValsNew(i) + (gamma/2)*edgePotentialsNew(i);
                 energyVals(i) = energyValNew;
                 totalEnergy = totalEnergy + energyValNew; 
             end
             
         end
        
        function [xVals,yVals] = calcNewXYVals(oldXVals, oldYVals,gradientEnergX, gradientEnergY,stepSize,radius,xCenter,yCenter)
            n = length(oldXVals);
            xVals = zeros(1,n);
            yVals = zeros(1,n);
            
            for i=1:n
                xVal = oldXVals(i) + gradientEnergX(i)*stepSize;
                yVal = oldYVals(i) + gradientEnergY(i)*stepSize;
                
                if snakeModel.isWithinCircle(radius,xCenter,yCenter,xVal,yVal)
                    xVals(i) = xVal;
                    yVals(i) = yVal;
                else
                    xVals(i) = oldXVals(i);
                    yVals(i) = oldYVals(i);
                end
                
            end
        end
                
        function withinCircle = isWithinCircle(radius,xCenter,yCenter,xVal,yVal)
            % Function to check if new values are within the init
            % circle/snake to prevent values outside of the image range
            dist = sqrt(( xVal - xCenter)^2 + (yVal - yCenter)^2);
            if dist < radius
                withinCircle = true;
            else
                withinCircle = false;
            end
            
        end
    end
    
    methods (Access= public)
        
         function [this,gradMag,newXVals,newYVals] = minimizeEnergy(this,stepSize)
             % function to minimize the snakes energy
             %  - calcs all needed derivatives and gradient Magnitude of
             %    edge image
             %  - calcs the gradient values
             %  - calcs the new x and y values for the snake contour
             %  - calcs new energy of snake
             
             %% Calc Gradient of EnergyFunc
             %2nd derivative of length/tension
             [secDerivX, secDerivY] = snakeModel.calcSecondDerivative(this.xVals,this.yVals);
             
             %4th derivative of stiffness
             [fourthDerivX, fourthDerivY] = snakeModel.calcFourthDerivative(this.xVals,this.yVals);
             
             %Gradient magnitude values of ImageEnergy (--> of edge image)
             [imageGradVals,~,gradMag] = snakeModel.calcGradientMagnitude(this.edgeImage, this.xVals, this.yVals);
             
             %new gradient values
             [gradientEnergyX,gradientEnergyY] = snakeModel.calcGradientVals(secDerivX,secDerivY,... 
                    fourthDerivX,fourthDerivY,imageGradVals, this.alpha, this.beta, this.gamma);  
                
             %new x,y coordinates
             [newXVals,newYVals] = snakeModel.calcNewXYVals(this.xVals,this.yVals, gradientEnergyX,gradientEnergyY,stepSize, this.radius, this.xCenter,this.yCenter);
             
             %new energy of snake
             [newEnergyVals,totalEnergyTmp] = snakeModel.calcEnergyVals(newXVals,newYVals,...
                 this.edgeImage,this.alpha,this.beta,this.gamma);
            
             %Set new values
             this.totalEnergy = totalEnergyTmp;
             this.energyVals = newEnergyVals;             
             this.xVals = newXVals;
             this.yVals = newYVals;  
             
         end
        
    end
end

