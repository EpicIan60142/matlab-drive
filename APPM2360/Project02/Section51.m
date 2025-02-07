%% Section 5.1: Image Translation
clc; clear; close all;


%% Problem 1

image = imread('rectangle.jpg');

[redMatrix, greenMatrix, blueMatrix, imageGray] = grayscale(image);

figure

hold on

sgtitle('Exracting Image Information and Grayscaling')

subplot(3,3,1)
colormap('gray')
imagesc(uint8(redMatrix));
subtitle('Red intensities')

subplot(3,3,3)
imagesc(uint8(greenMatrix));
subtitle('Green intensities')

subplot(3,3,7)
imagesc(uint8(blueMatrix));
subtitle('Blue intensities')

subplot(3,3,9)
imagesc(uint8(imageGray));
subtitle('Grayscale image')

subplot(3,3,5)
imagesc(image)
subtitle('Original image')

imwrite(uint8(imageGray),"grayRectangle.jpg");
imwrite(uint8(redMatrix),"redRectangle.jpg");
imwrite(uint8(greenMatrix),"greenRectangle.jpg");
imwrite(uint8(blueMatrix),"blueRectangle.jpg");

%% Problem 2

figure
colormap('gray')

sgtitle('Increased Grayscale Exposure of an Image')

imageWhitedOut = imageGray + 100;

subplot(1,2,2)
imagesc(uint8(imageWhitedOut));
subtitle('Increased exposure');

subplot(1,2,1)
imagesc(uint8(imageGray));
subtitle('Original grayscale')

imwrite(uint8(imageWhitedOut),"increasedExposureRectangle.jpg");

%% Problem 3

figure

sgtitle('Modifying Color Intensities of an Image')

redMatrix2 = 0*redMatrix;
greenMatrix2 = greenMatrix;
blueMatrix2 = blueMatrix + 80;

image2(:,:,1) = redMatrix2;
image2(:,:,2) = greenMatrix2;
image2(:,:,3) = blueMatrix2;

subplot(1,2,2)
imagesc(uint8(image2))
subtitle('Modified')

subplot(1,2,1)
imagesc(uint8(image))
subtitle('Original')

imwrite(uint8(image2),"modifiedRectangle.jpg");

%% Problem 4

A = [1 2; 3 4]
E = [0 1; 1 0];
A*E
E*A

% E is the matrix 0 1, and it should be multiplied to the right of A.
%                 1 0

%% Problem 5

figure

sgtitle('Horizontal and Vertical Image Translation')

vertShiftImage = vertTranslate(image, 306);
horizShiftImage = horizTranslate(image, 306);

subplot(1,3,3)
imagesc(uint8(vertShiftImage))
subtitle('Vertical shift')

subplot(1,3,1)
imagesc(uint8(image))
subtitle('Original')

subplot(1,3,2)
imagesc(uint8(horizShiftImage))
subtitle('Horizontal shift')

imwrite(uint8(vertShiftImage),"vertRectangle.jpg");
imwrite(uint8(horizShiftImage),"horizRectangle.jpg");

%% Problem 6

figure

sgtitle('Combined Translation of an Image')

[vertShiftImage2, vertMatrix] = vertTranslate(image, 230);
[horizShiftImage2, horizMatrix] = horizTranslate(vertShiftImage2, 306);

subplot(2,3,1)
imagesc(uint8(image))
subtitle('Original')

subplot(2,3,2)
imagesc(uint8(vertShiftImage2))
subtitle('Vertical shift')

subplot(2,3,3)
imagesc(uint8(horizShiftImage2))
subtitle('Vertical and horizontal shift')

subplot(2,3,4)
spy(vertMatrix)
subtitle('Vertical transformation matrix')

subplot(2,3,5)
spy(horizMatrix)
subtitle('Horizontal transformation matrix')

imwrite(uint8(vertShiftImage2),"vertRectangle2.jpg");
imwrite(uint8(horizShiftImage2),"combinedShiftRectangle.jpg");

%% Problem 7

figure

sgtitle('Flipping an Image')

[flippedImage, flipMatrix] = flipImage(image);

subplot(1,3,1)
imagesc(uint8(image));
subtitle('Original')

subplot(1,3,2)
imagesc(uint8(flippedImage));
subtitle('Flipped image')

subplot(1,3,3)
spy(flipMatrix)
subtitle('Flip transformation')

imwrite(uint8(flippedImage),"flippedRectangle.jpg");

%% Problem 8

%Transposing an image should flip it along the main diagonal of the matrix.

figure
colormap('gray')

sgtitle('Transposing an Image')

transposedImage(:,:,1) = redMatrix';
transposedImage(:,:,2) = greenMatrix';
transposedImage(:,:,3) = blueMatrix';

subplot(1,3,1)
imagesc(image);
subtitle('Original');

subplot(1,3,2)
imagesc(transposedImage)
subtitle('Color transposed')

subplot(1,3,3)
imagesc(imageGray')
subtitle('Grayscale transposed')

%Color doesn't look as expected, grayscale does.

imwrite(uint8(transposedImage),"colorTransposedRectangle.jpg");
imwrite(uint8(imageGray'),"grayscaleTransposedRectangle.jpg");

%% Problem 9

figure

sgtitle('Cropping an Image')

[croppedImage, E1, E2, border] = crop(image, 50);

subplot(2,2,1)
imagesc(uint8(image));
subtitle('Original')

subplot(2,2,2)
imagesc(uint8(croppedImage))
subtitle(['Cropped image, border = ',num2str(border),' px'])

subplot(2,2,3)
spy(E1)
subtitle('Top and bottom bands')

subplot(2,2,4)
spy(E2)
subtitle('Left and right bands')

imwrite(uint8(croppedImage),"croppedRectangle.jpg");
