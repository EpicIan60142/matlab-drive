function [redMatrix, greenMatrix, blueMatrix, grayImage] = grayscale(image)
% Takes in an image and converts it into a double grayscaled image
%   Inputs: 8-bit integer image
%   Outputs: double red intensity matrix, double green intensity matrix, double blue intensity
%   matrix, double grayscale image
    redMatrix = double(image(:,:,1));
    greenMatrix = double(image(:,:,2));
    blueMatrix = double(image(:,:,3));
    
    imageDouble = double(image);
    grayImage = imageDouble(:,:,1)/3.0 + imageDouble(:,:,2)/3.0 + imageDouble(:,:,3)/3.0;
end

