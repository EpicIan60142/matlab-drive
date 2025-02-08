function [translatedImage, T] = vertTranslate(image, shift)
% Shifts an image vertically some amount
%   Inputs: Image to shift, shift size
%   Outputs: Shifted image

doubleImage = double(image);

for k = 1:3
    [m,n] = size(doubleImage(:,:,k));
    r = shift;  
    E = eye(m);
    T = zeros(m);
    % fill in the first r rows of T with the last r rows of E
    T(1:r,:) = E(m-(r-1):m,:);
    % fill in the rest of T with the first part of E
    T(r+1:m,:) = E(1:m-r,:);
    translatedImage(:,:,k) = T*doubleImage(:,:,k);
end

end

