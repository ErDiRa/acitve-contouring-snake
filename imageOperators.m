classdef imageOperators
    
    methods(Static)
        
        function res = gaussianFilter(imageData)
            sig = 1; 
            siz = 4*sig + 1;
            hG = fspecial('gaussian', siz, sig);
            res = uint8(conv2(imageData, hG, 'same'));
        end
        
        function res = anisotropicFilter(imageData)
            imageData = double(imageData);
            res = imdiffusefilt(imageData);
        end

        function res = cannyFilter(imageData)
             res = edge(imageData,'canny');
        end

        function res = detectEdgesLaplace(imageData)
            tau = 4;
            laPlaceFilter = [0 -1 0; -1 4 -1; 0 -1 0];
            laPlaceFilter = tau.*laPlaceFilter;
            firstStep = conv2(imageData,laPlaceFilter);
            res = uint8(firstStep);
        end

        function res = detectEdgesXY(imageData)
            hx = [-1,0,1];
            hy = hx';
            filterMask = hy*hx;
            edgeImage = conv2(imageData,filterMask);
            res = uint8(edgeImage);
        end

        function res = detectEdgesSobel(imageData)
            % filter is blurring and detecting edges
            sobelX = [-1 0 1; -2 0 2; -1 0 1];
            sobelY = [-1 -2,-1; 0 0 0 ; 1 2 1];
            convX = conv2(imageData, sobelX);
            convY = conv2(imageData, sobelY);
            gradient = sqrt((convX).^2 + (convY).^2);
            res = uint8(gradient);
        end
        
        function output = convertToGrey(image)
            alpha= 1/3;beta = 1/3;gamma = 1/3;
            [height,width,depth] = size(image);
            output = zeros(height,width);
            output = alpha * image(:,:,1) +  beta * image(:,:,2) + gamma * image(:,:,3);
        end

        function output = performThresholding(image, threshold)
            [height, width] = size(image);
            output = zeros(height,width);
            for i=1:height
                for j = 1: width
                    if image(i,j) > threshold
                      output(i,j) = 255;
                    else
                        output(i,j) = 0;
                    end
                end
            end

        end
        
        function image_stretched = performWindowing(image,gmin,gmax,gmax_win)
            % g < gmin --> f(g) = 0
            % gmin < g < gmax -- f(g) = g'max * (g-gmin)/(gmax-gmin)
            % g > gmax = 255

            image = double(image);
            [M,N] = size(image);
            image_stretched = zeros(M,N);

            for i=1:M
                for j=1:N
                    greyVal = image(i,j);
                    if greyVal < gmin
                        newGreyVal = 0;
                    elseif ((gmin <= greyVal) && (greyVal<= gmax))
                        newGreyVal = round(gmax_win * ((greyVal -gmin)/(gmax-gmin)));
                    elseif greyVal > gmax
                        newGreyVal = gmax_win;
                    end
                    image_stretched(i,j) = newGreyVal;
                end

            end
    
            image_stretched = uint8(image_stretched);

        end

        function image_sharpened = performSharpening(input_image,factor)

            blurred = imgaussfilt(input_image);
            sharpness = input_image - blurred;

            image_sharpened = input_image + factor * sharpness;

        end
        
        function convertedGrayVals= convertGreyValsToInt8(greyValuesDicom)
            %CONVERTGREYVALS scale dicom image values to 0-255
            greyValuesDicom = double(greyValuesDicom);
            greyValuesDicom = greyValuesDicom / max(greyValuesDicom(:));
            convertedGrayVals= uint8(255* greyValuesDicom);
        end
        
        function image_medFiltered = medianFilter(imageData)
            image_medFiltered = medfilt2(imageData);
        end
    end
end

