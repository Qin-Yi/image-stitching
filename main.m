%% Clear all
clc; close all; clc;

%% Load a list of images
imgList = dir('./data/tree*.jpg');
saveFileName = 'tree.jpg';

IMAGES = cell(1, length(imgList));
for i = 1 : length(imgList)
    IMAGES{i} = imread(['./data/' imgList(i).name]);
    %% Resize to make memory efficient
    if max(size(IMAGES{i})) > 1000 || length(imgList) > 10
        IMAGES{i} = imresize(IMAGES{i}, 0.2); % 0.1
    end
end

%% Feature detection
DESCRIPTOR = cell(1, length(imgList));
POINT_IN_IMG = cell(1, length(imgList));
for i = 1 : length(imgList)
    [feature,imp] = detect_features(IMAGES{i});
    POINT_IN_IMG{i} = feature(:, 1:2);
    pointInPyramid = feature(:, 4:5);
    DESCRIPTOR{i} = SIFTDescriptor(imp, pointInPyramid, feature(:,3));
end

%% Compute Transformation
TRANSFORM = cell(1, length(imgList)-1);
for i = 1 : (length(imgList)-1)
    TRANSFORM{i} = RANSACFit(POINT_IN_IMG{i}, POINT_IN_IMG{i+1}, DESCRIPTOR{i}, DESCRIPTOR{i+1});
end

%% Make Panoramic image
MultipleStitch(IMAGES, TRANSFORM, saveFileName);
imshow(imread(saveFileName));


%% scale space feature detector based upon difference of gaussian filters.
% selects features based upon their maximum response in scale space
function [features,imp]  = detect_features(img)
%   recommended parameter values from Stanford CS131 
    scale = 1.5;
    threshold = 3;
    radius = 4;

    if size(img,3) > 1
        img = rgb2gray(img);
    end
    
    % Computation of the maximum number of levels:
    Lmax = floor(min(log(2*size(img)/12)/log(scale)));
    
    % Build image pyramid and difference of gaussians filter response pyramid
    [pyr,imp] = build_pyramid(img,Lmax,scale);  
    
    % Get the feature points
    pts = find_features(pyr,scale,threshold,radius);

    % classify points and create sub-pixel and sub-scale adjustments 
    [features] = refine_features(img,pyr,scale,imp,pts,radius);
end

%% build scaled image pyramid and difference of gaussians pyramid
function [pyr,imp] = build_pyramid(img,levels,scale)
    % expand to retain spatial frequencies
    img2 = resample_bilinear(img,1/2);         

    % variance for laplacian filter
    sigma=1.5;    
    % variance for downsampling
    sigma2=1.5;  
    sig_delta = 1/levels;

    for i=1:levels        
         if i==1
             img3 = img2;
             % slightly filter bottom level
             img2 = imgaussfilt(img2,.5); 
         end

         imp{i}=img2;
         % calculate difference of gaussians
         A = imgaussfilt(img2,sigma);   
         B = imgaussfilt(A,sigma);
         pyr{i} = A-B;

         if i==1
             img2 = img3;
         else
             % anti-alias for downsampling
             B = imgaussfilt(img2,sigma2);    
             B = imgaussfilt(B,sigma2);
         end

         img2 = resample_bilinear(B,scale);           
    end
end

%% resamples a 2D matrix by the ratio using bilinear interpolation
% the 1,1 entry of the matrix is always duplicated. 
function img2 = resample_bilinear(img, ratio)
    img=double(img);
    [h,w]=size(img); 
    
    % create vectors of X and Y values for new image
    [y,x]=meshgrid( 1:ratio:h-1, 1:ratio:w-1 );  
    [h2,w2] = size(x);                 
    % convert to vectors
    x = x(:); y = y(:);
    %calculate alphas and betas for each point
    alpha = x - floor(x);   
    beta = y - floor(y);
    fx = floor(x); fy = floor(y);

    % create index for neighboring pixels
    inw = fy + (fx-1)*h;    
    ine = fy + fx*h;
    isw = fy+1 + (fx-1)*h;
    ise = fy+1 + fx*h;
    
    %interpolate
    img2 = (1-alpha).*(1-beta).*img(inw) + (1-alpha).*beta.*img(isw) + ...
           alpha.*(1-beta).*img(ine) + alpha.*beta.*img(ise);

    img2 = reshape(img2,h2,w2);                 
    img2 = img2';
end


