classdef snakeModel
    
    properties
        alpha
        beta
        gamma
        finalXVals
        finalYVals
        xVals
        yVals
        energyVals
        totalEnergy
        energyValsInit
        totalEnergyInit
        edgeImage
        
    end
    
    methods(Static)
        
        function snake = create(alpha,beta,gamma,xVals,yVals,imageData)
            snake = snakeModel;
            snake.alpha = alpha;
            snake.beta = beta;
            snake.gamma = gamma;
            snake.xVals = xVals;
            snake.yVals = yVals;
            snake.edgeImage = snakeModel.prepareEdgeImage(imageData);
            [snake.energyValsInit, snake.totalEnergyInit ] = ...
                snakeModel.calcInitEnergyVals(xVals,yVals,snake.edgeImage, snake.alpha,snake.beta,snake.gamma);
            snake.totalEnergy = snake.totalEnergyInit;
            snake.energyVals = snake.energyValsInit;
        end
        
    end
    
    methods(Access = private, Static)
        
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
% 
             canny_Image = imageOperators.cannyFilter(blurred_Image);
             canny_Image_blurred = imageOperators.anisotropicFilter(canny_Image);
             subplot(2,3,5)
             imshow(canny_Image)
             title('canny filtered')

            %% Inverting
            %edgeImage = imcomplement(edge_Image);%if snake hits an edge potential gets zero
            edgeImage = canny_Image_blurred;
            subplot(2,3,4)
            imshow(edgeImage)
            title('edge Image')
            %TODO: the edge image are forces maybe try to visualize that
         end 
        
        function [gradMagValues,gradMagDirection,gradMag] = calcGradientMagnitude(edgeImage,xVals,yVals)
            [gradMag, gradMagDir] = imgradient(edgeImage);
            gradMagValues = zeros(1,length(xVals));
            gradMagDirection = zeros(1, length(yVals));
            for i=1:length(xVals)
                xpos = round(xVals(i));
                ypos = round(yVals(i));
                gradMagValues(i) = gradMag(ypos,xpos);
                gradMagDirection = gradMagDir(ypos,xpos);
            end
            
            gradMag = uint8(gradMag);
            gradMagValues = uint8(gradMagValues);
            
        end
        
        function [secDerivX, secDerivY] = calcSecondDerivative(xVals,yVals)
           
             n = length(xVals);
             secDerivX = zeros(1,n);
             secDerivY = zeros(1,n);
             for i = 1:n
                 
                 if i == 1
                    v_prev = [xVals(n-1); yVals(n-1)]; %because at index n is the same value as 1 (close circle)
                    v = [xVals(i); yVals(i)];
                    v_next = [xVals(i+1);yVals(i+1)];

                    %2nd derivative
                    secDerivX(i) = v_prev(1) - 2*v(1) + v_next(1); %%times stepsize?times stepsize?
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
             n = length(xVals);
             fourthDerivX = zeros(1,n);
             fourthDerivY = zeros(1,n);
             for i = 1:n
                 
                 if i == 1
                    v_pprev = [xVals(n-2); yVals(n-2)];
                    v_prev = [xVals(n-1); yVals(n-1)]; %because at index n is the same value as 1 (close circle)
                    v = [xVals(i); yVals(i)];
                    v_next = [xVals(i+1);yVals(i+1)];
                    v_nnext = [xVals(i+2);yVals(i+2)];

                    %4th derivative
                    fourthDerivX(i) = v_pprev(1) - 4* v_prev(1) + 6 * v(1) - 4* v_next(1) + v_nnext(1);
                    fourthDerivY(i) = v_pprev(2) - 4* v_prev(2) + 6 * v(2) - 4* v_next(2) + v_nnext(2);
                    
                 elseif i == 2
                    v_pprev = [xVals(n-1); yVals(n-1)];
                    v_prev = [xVals(1); yVals(1)]; %because at index n is the same value as 1 (close circle)
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
                elseif i ==n
                    tensionValues(i) = tensionValues(1);
                end
                
            end
            
            totalTension = sum(tensionValues);
         
         end
         
        function [stiffnessVals,stiffnessTotal] = calcStiffness(xVals,yVals)
             n = length(xVals);
             stiffnessVals = zeros(1, n);
             
             for i = 1:n
                 
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
            %%return image values of edge detected image
            imageEnergyVals = zeros(1, length(xVals));
            for i=1:length(xVals)
                xpos = round(xVals(i));
                ypos = round(yVals(i));
                
                imageEnergyVals(i) = double((edgeImage(ypos,xpos)))^2;
            end
            
            imageEnergyTotal = sum(imageEnergyVals);
        end
        
        function [gradientEnergyX,gradientEnergyY] = calcGradientEnergies(secDerivX,secDerivY, fourthDerivX,...
                fourthDerivY,imageGradVals,alpha,beta,gamma)
            
            n = length(secDerivX);
            gradientEnergyX = zeros(1,n);
            gradientEnergyY = zeros(1,n);
            
            for i=1:n
                gradientEnergyX(i) = alpha*secDerivX(i) + beta*fourthDerivX(i) + gamma*double(imageGradVals(i));
                gradientEnergyY(i) = alpha*secDerivY(i) + beta*fourthDerivY(i) + gamma*double(imageGradVals(i));
            end
            
            
        end
        
        function [energyVals,totalEnergy] = calcInitEnergyVals(xVals,yVals,edgeImage,alpha,beta,gamma)
             
             [tensionValsNew, ~] = snakeModel.calcTension(xVals,yVals);
  
             n = length(tensionValsNew(1,:));
             
             [stiffnessValsNew, ~] = snakeModel.calcStiffness(xVals, yVals);
             
             [edgePotentialsNew,~] = snakeModel.calcImageForces(xVals,yVals,edgeImage);
             
             energyVals = zeros(1,n);
             totalEnergy = 0;
             for i=1:n
                 energyVals(i) = (alpha/2)*tensionValsNew(i) + (beta/2)*stiffnessValsNew(i) + (gamma/2)*edgePotentialsNew(i);
             end
             
         end
        
        function [newXVals,newYVals,energyVals,totalEnergy] = calcEnergyVals(xVals,yVals,oldXVals,oldYVals,...
                totalEnergySnake,oldEnergyVals,edgeImage,alpha,beta,gamma)
             
             [tensionValsNew, ~] = snakeModel.calcTension(xVals,yVals);
             [tensionValsOld, ~] = snakeModel.calcTension(oldXVals,oldYVals);
             n = length(tensionValsNew(1,:));
             newXVals = zeros(1,n);
             newYVals = zeros(1,n);
             
             [stiffnessValsNew, ~] = snakeModel.calcStiffness(xVals, yVals);
             [stiffnessValsOld, ~] = snakeModel.calcStiffness(oldXVals, oldYVals);
             
             [edgePotentialsNew,~] = snakeModel.calcImageForces(xVals,yVals,edgeImage);
             [edgePotentialsOld,~] = snakeModel.calcImageForces(oldXVals,oldYVals,edgeImage);
             
             energyVals = oldEnergyVals;
             totalEnergy = 0;
             for i=1:n
                 energyValNew = (alpha/2)*tensionValsNew(i) + (beta/2)*stiffnessValsNew(i) + (gamma/2)*edgePotentialsNew(i);
                 energyValOld = (alpha/2)*tensionValsOld(i) + (beta/2)*stiffnessValsOld(i) + (gamma/2)*edgePotentialsOld(i);
                 
                 if (energyValNew < energyValOld) 
                     energyVals(i) = energyValNew;
                     if sum(energyVals) < sum(oldEnergyVals) % remove this, because when value is smaller its sum is also smaller
                         newXVals(i) = xVals(i);
                         newYVals(i) = yVals(i);
                     else
                         energyVals(i) = energyValOld;
                         newXVals(i) = oldXVals(i);
                         newYVals(i) = oldYVals(i);
                     end
                      totalEnergy = totalEnergy + energyValNew;
                 else
                     energyVals(i) = energyValOld;
                     newXVals(i) = oldXVals(i);
                     newYVals(i) = oldYVals(i);
                     totalEnergy = totalEnergy + energyValOld;
                 end
             end
             
         end
        
        function [xVals,yVals] = calcNewXYVals(oldXVals, oldYVals,gradientEnergX, gradientEnergY,stepSize)
            n = length(oldXVals);
            xVals = zeros(1,n);
            yVals = zeros(1,n);
            
            for i=1:n
                %check if new xVal is within range otherwise take the old
                %one
                xVals(i) = oldXVals(i) + gradientEnergX(i)*stepSize;
                yVals(i) = oldYVals(i) + gradientEnergY(i)*stepSize;
% %                 startIdx = i; endIdx = n-i; if i < (n-1)/2
% %                     startIdx = i; endIdx = n-i; [xInRange,yInRange] =
% %                     snakeModel.withinRange(xVal,yVal,
% %                     initXVals,initYVals,startIdx,endIdx,radius);
% %                 else
% %                     startIdx = n-i; endIdx = ; [xInRange,yInRange] =
% %                     snakeModel.withinRange(xVal,yVal,
% %                     initXVals,initYVals,endIdx,startIdx,radius);
% %                 end if xInRange
% %                     xVals(i) = xVal;
% %                 else
% %                     xVals(i) = oldXVals(i);
% %                 end
% %                 
% %                 if yInRange
% %                     yVals(i) = yVal;
% %                 else
% %                     yVals(i) = oldYVals(i);
% %                 end
                
            end
        end
       
        
    end
    
   
    
    methods (Access= public)
        
         function [this,gradMag,newXVals,newYVals] = minimizeEnergy(this,stepSize)
             % therefore we need to derive the energy function
             % https://en.wikipedia.org/wiki/Finite_difference use finite
             % difference to get the derivatives
             % necessary: 
             %           2nd derivative of the length
             %           4th derivative of the stiffness
             %           Gradient magnitude values of the Image Energy
             %           (edge image)
             
             
        
             %% Calc Gradient of EnergyFunc
             %2nd derivative of length
             [secDerivX, secDerivY] = snakeModel.calcSecondDerivative(this.xVals,this.yVals);
             
             %4th derivative of stiffness
             [fourthDerivX, fourthDerivY] = snakeModel.calcFourthDerivative(this.xVals,this.yVals);
             
             %Gradient magnitude values of ImageEnergy (--> of edge image)
             [imageGradVals,~,gradMag] = snakeModel.calcGradientMagnitude(this.edgeImage, this.xVals, this.yVals);
             
             [gradientEnergyX,gradientEnergyY] = snakeModel.calcGradientEnergies(secDerivX,secDerivY,... 
                    fourthDerivX,fourthDerivY,imageGradVals, this.alpha, this.beta, this.gamma);  
             
             [newXVals,newYVals] = snakeModel.calcNewXYVals(this.xVals,this.yVals, gradientEnergyX,gradientEnergyY,stepSize);
             [newXVals,newYVals,newEnergyVals,totalEnergyTmp] = snakeModel.calcEnergyVals(newXVals,newYVals,this.xVals,this.yVals,...
                 this.totalEnergy,this.energyVals,this.edgeImage,this.alpha,this.beta,this.gamma);
            
             this.totalEnergy = totalEnergyTmp;
             this.energyVals = newEnergyVals;
             if totalEnergyTmp < this.totalEnergy
                 this.totalEnergy = totalEnergyTmp;
                 this.finalXVals = newXVals;
                 this.finalYVals = newYVals;
             end
             
             this.xVals = newXVals;
             this.yVals = newYVals;  
             
         end
        
    end
end

