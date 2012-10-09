function im = give_im(imPath,invert,imMean)

% Load images
im = imread(imPath);

% Adjust grayscale values
im   = (imadjust(im));

%img = adapthisteq(img,'clipLimit',0.02,'Distribution','rayleigh');

% Subtract background
if nargin > 2
    warning off
    im = imsubtract(imadjust(imMean),im);
    warning on
    
    if ~invert
        im = imcomplement(im);
    end
else
    if invert
        im = imcomplement(im);
    end
end

% Enhance contrast
im = adapthisteq(im);