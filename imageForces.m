function [potVal,image] = imageForces(imageData)
%IMAGEFORCES draws Potential towards an edge of image
%   edge detection: x -> x + ddI(x) + dI(x)

%blurred_Image = gaussianFilter(imageData);
blurred_Image = anisotropicFilter(imageData);
edge_Image = detectEdges(blurred_Image);
image = imcomplement(edge_Image); %if snake hits an edge potential gets zero
potVal = norm(double(image)); %I think this must be the potvalue of the snake and not the image

end

function res = gaussianFilter(imageData)
    sig = 1; siz = 6*sig + 1;
    hG = fspecial('gaussian', siz, sig);
    res = conv2(imageData, hG, 'same');
end

function res = anisotropicFilter(imageData)
     res = imdiffusefilt(imageData);
end


function res = detectEdges(imageData)
    tau = 2;
    
    laPlaceFilter = [0 -1 0; -1 4 -1; 0 -1 0];
    laPlaceFilter = tau.*laPlaceFilter;
    firstStep = conv2(imageData,laPlaceFilter);
    
    hx = [-1,0,1];
    hy = hx';
    filterMask = hy*hx;
    secondStep = conv2(firstStep,filterMask);
    
    res = uint8(secondStep);
    
end
