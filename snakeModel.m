classdef snakeModel
    
    properties
        alpha
        beta
        gamma
        xVals
        yVals
        energyVal
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
        end
        
    end
    
    methods(Access = private,Static)
        
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

%             canny_Image = imageOperators.cannyFilter(blurred_Image);
%             subplot(2,3,5)
%             imshow(canny_Image)
%             title('canny filtered')

            %% Inverting
            %this.edgeImage = imcomplement(edge_Image);%if snake hits an edge potential gets zero
            edgeImage = edge_Image;
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
                gradMagValues(i) = gradMag(xVals(i),yVals(i));
                gradMagDirection = gradMagDir(xVals(i),yVals(i));
            end
            
        end
        
        function [secDerivX, secDerivY] = calcSecondDerivative(xVals,yVals)
             %ToDO do not forget h
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
                 elseif i == n
                      secDerivX(i) = secDerivX(1);
                      secDerivY(i) = secDerivY(1);
                 elseif i < n 
                    
                    v_prev = [xVals(i-1); yVals(i-1)];
                    v = [xVals(i); yVals(i)];
                    v_next = [xVals(i+1);yVals(i+1)];

                    %2nd derivative
                    secDerivX(i) = v_prev(1) - 2*v(1) + v_next(1); %times stepsize?
                    secDerivY(i) = v_prev(2) - 2*v(2) + v_next(2);%times stepsize?
                 end
   
             end
         end
         
        function [fourthDerivX, fourthDerivY] = calcFourthDerivative(xVals,yVals)
             %ToDO do not forget h
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
                    fourthDerivX(i) = v_pprev(1) - 4* v_prev(1) + 6 * v(1) - 4* v_next(1) + v_nnext(1);
                    fourthDerivY(i) = v_pprev(2) - 4* v_prev(2) + 6 * v(2) - 4* v_next(2) + v_nnext(2);     
                    
                 elseif i == n-1 
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
                 else
                    v_pprev = [xVals(i-2); yVals(i-2)];
                    v_prev = [xVals(i-1); yVals(i-1)];
                    v = [xVals(i); yVals(i)];
                    v_next = [xVals(i+1);yVals(i+1)];
                    v_nnext = [xVals(i+2);yVals(i+2)];
                    
                    %4th derivative
                    fourthDerivX(i) = v_pprev(1) - 4* v_prev(1) + 6 * v(1) - 4* v_next(1) + v_nnext(1);
                    fourthDerivY(i) = v_pprev(2) - 4* v_prev(2) + 6 * v(2) - 4* v_next(2) + v_nnext(2);
                 end
                 
                
   
             end
        end
         
        function [tensionValues,totalTension] = calcTension(xVals,yVals)
            
            n = length(xVals);
            tensionValues = zeros(1,n);
         
                
            for i=1:n
                
                if i < n
                
                    v = [xVals(i); yVals(i)];
                    v_next = [xVals(i+1);yVals(i+1)];

                    %1st derivative (forward difference)
                    deriveX = v_next(1) - v(1);
                    deriveY = v_next(2) - v(2);


                    tensionValues(i) = (sqrt(deriveX^2 + deriveY^2))^2;
               
                else
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
                imageEnergyVals(i) = double((edgeImage(xVals(i),yVals(i))))^2; %TODO check out why it is y,x not x,y
            end
            
            imageEnergyTotal = sum(imageEnergyVals);
        end
        
    end
    
    methods
        
         function [tensionVals, stiffnessVals, edgePotentials, energyVals, totalEnergy] = calcEnergyVals(this)
             
             
             [tensionVals, tensionTotal] = snakeModel.calcTension(this.xVals, this.yVals);
             n = length(tensionVals(1,:));
             
             [stiffnessVals, stiffnessTotal] = snakeModel.calcStiffness(this.xVals, this.yVals);
             
             [edgePotentials,edgePotentialTotal] = snakeModel.calcImageForces(this.xVals,this.yVals,this.edgeImage);
             
             energyVals = zeros(1,n);
    
             for i=1:n
                 energyVals(i) = (this.alpha/2)*tensionVals(1,i) + (this.beta/2)*stiffnessVals(1,i) + (this.gamma/2)*edgePotentials(1,i);
             end
            
             totalEnergy = (this.alpha/2)*tensionTotal + (this.beta/2)*stiffnessTotal + (this.gamma/2)*edgePotentialTotal;
             
         end
         
         function [secDerivX, secDerivY, fourthDerivX, fourthDerivY,imageGradVals,gradMag] = minimizeEnergy(this,stepSize)
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
             
         end
        
    end
end

