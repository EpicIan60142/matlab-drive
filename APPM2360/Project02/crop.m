function [croppedImage, E1, E2, border] = crop(image, b)
% Crops a specified border size around an image
%   Inputs: An image, border size
%   Outputs: Cropped image, transformation matrix for top and bottom cropping,
%   transformation matrix for left and right cropping

redMatrix = double(image(:,:,1));
greenMatrix = double(image(:,:,2));
blueMatrix = double(image(:,:,3));

[m,n] = size(redMatrix);
border = b;

E1 = eye(m);
E2 = eye(n);

%Top and bottom bands
E1(1:b,:) = 0;
E1(m-b:m,:) = 0;

%Left and right bands
E2(:,1:b) = 0;
E2(:,n-b:n) = 0;

croppedImage(:,:,1) = E1*redMatrix*E2;
croppedImage(:,:,2) = E1*greenMatrix*E2;
croppedImage(:,:,3) = E1*blueMatrix*E2;


end

