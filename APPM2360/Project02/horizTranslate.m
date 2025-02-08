function [translatedImage, T] = horizTranslate(image, shift)
% Shifts an image horizontally some amount
%   Inputs: Image to shift, shift size
%   Outputs: Shifted image

doubleImage = double(image);

for k = 1:3
    [m,n] = size(doubleImage(:,:,k));
    r = shift;
    E = eye(n);
    T = zeros(n);
    % fill in the first r rows of T with the last r rows of E
    T(:,1:r) = E(:,n-(r-1):n);
    % fill in the rest of T with the first part of E
    T(:,r+1:n) = E(:,1:n-r);
    translatedImage(:,:,k) = doubleImage(:,:,k)*T;
end
    
end
