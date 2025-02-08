function [flippedImage, T] = flipImage(image)
% Flips an image across the horizontal axis
%   Inputs: An image
%   Output: Flipped image, image transformation matrix

doubleImage = double(image);

for k = 1:3
    [m,n] = size(doubleImage(:,:,k));
    E = eye(m);
    T = zeros(m);

    T(:,:) = E(end:-1:1,:); %Assign T values corresponding to indexing 
                            %backwards through the identity matrix, 
                            %effectively flipping it
    
    flippedImage(:,:,k) = T*doubleImage(:,:,k);
    
end
    
end

