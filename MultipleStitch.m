function Pano = MultipleStitch( IMAGES, TRANS, fileName )
%   This function stitches multiple images together and outputs the panoramic stitched image
%   with a chain of input images and its corresponding transformations. 
%   
%   Given a chain of images:
%       I1 -> I2 -> I3 -> ... -> Im
%
%   We choose the middle image as the reference image, and the outputed
%   panorama is in the same coordinate system as the reference image.

if length(IMAGES) ~= length(TRANS)+1
    error('Number of images does not match the number of transformations.');
end

%% Outbounds of panorama image
outBounds = zeros(2,2);
outBounds(1,:) = Inf;
outBounds(2,:) = -Inf;

%% Choose reference image Iref
ref_idx = ceil(median(1:length(IMAGES)));

%% Estimate the largest possible panorama size
[nrows, ncols, ~] = size(IMAGES{1});
nrows = length(IMAGES) * nrows;
ncols = length(IMAGES) * ncols;

% projective transformation from IMAGES{i} to the reference image
imageToRefTrans = cell(1, length(IMAGES));

% Initialize imageToRefTrans to contain the identity transform.
for idx = 1:length(imageToRefTrans)
    imageToRefTrans{idx} = eye(3);
end

%% Find the correct transformations used for images on the left side of Iref
for idx = ref_idx-1:-1:1
    imageToRefTrans{idx} = makeTransformToReferenceFrame(TRANS, idx, ref_idx);
    T = imageToRefTrans{idx};
    tmpBounds = findbounds(maketform('projective', T'), [1 1; ncols nrows]);
    outBounds(1,:) = min(outBounds(1,:),tmpBounds(1,:));
    outBounds(2,:) = max(outBounds(2,:),tmpBounds(2,:));
end

%% Find the correct transformations used for images on the right side of Iref

for idx = ref_idx + 1 : length(imageToRefTrans)  
    imageToRefTrans{idx} = makeTransformToReferenceFrame(TRANS, idx, ref_idx);
    T = imageToRefTrans{idx};
    tmpBounds = findbounds(maketform('projective', T'), [1 1; ncols nrows]);
    outBounds(1,:) = min(outBounds(1,:),tmpBounds(1,:));
    outBounds(2,:) = max(outBounds(2,:),tmpBounds(2,:));
end

%% Stitch the Iref image.
XdataLimit = round(outBounds(:,1)');
YdataLimit = round(outBounds(:,2)');
Pano = imtransform( im2double(IMAGES{ref_idx}), maketform('projective', eye(3)), 'bilinear', ...
                    'XData', XdataLimit, 'YData', YdataLimit, ...
                    'FillValues', NaN, 'XYScale',1);
                
%% Transform the images from the left side of Iref using the correct transformations you computed
for idx = ref_idx-1:-1:1
    T = imageToRefTrans{idx};
    Tform = maketform('projective', T');
    AddOn = imtransform(im2double(IMAGES{idx}), Tform, 'bilinear', ...
                        'XData', XdataLimit, 'YData', YdataLimit, ...
                        'FillValues', NaN, 'XYScale',1);
    result_mask = ~isnan(Pano(:,:,1));
    temp_mask = ~isnan(AddOn(:,:,1));
    add_mask = temp_mask & (~result_mask);

    for c = 1 : size(Pano,3)
        cur_im = Pano(:,:,c);
        temp_im = AddOn(:,:,c);
        cur_im(add_mask) = temp_im(add_mask);
        Pano(:,:,c) = cur_im;
    end
end

%% Transform the images from the right side of Iref using the correct transformations you computed
for idx = ref_idx + 1 : length(imageToRefTrans)
    T = imageToRefTrans{idx};
    Tform = maketform('projective', T');
    AddOn = imtransform(im2double(IMAGES{idx}), Tform, 'bilinear', ...
                        'XData', XdataLimit, 'YData', YdataLimit, ...
                        'FillValues', NaN, 'XYScale',1);
    result_mask = ~isnan(Pano(:,:,1));
    temp_mask = ~isnan(AddOn(:,:,1));
    add_mask = temp_mask & (~result_mask);

    for c = 1 : size(Pano,3)
        cur_im = Pano(:,:,c);
        temp_im = AddOn(:,:,c);
        cur_im(add_mask) = temp_im(add_mask);
        Pano(:,:,c) = cur_im;
    end
end

%% Cropping the final panorama to leave out black spaces.
[I, J] = ind2sub([size(Pano, 1), size(Pano, 2)], find(~isnan(Pano(:, :, 1))));
upper = max(min(I)-1, 1);
lower = min(max(I)+1, size(Pano, 1));
left = max(min(J)-1, 1);
right = min(max(J)+1, size(Pano, 2));
Pano = Pano(upper:lower, left:right,:);

imwrite(Pano, fileName);

end

function T = makeTransformToReferenceFrame(i_To_iPlusOne_Transform, idx, ref_idx)
% T: A 3x3 homogeneous transformation matrix that 
% convert a point in the current frame into the corresponding point in the reference frame. 
% If the current frame and reference frame are not adjacent, T will need to be calculated.

if idx<ref_idx
    T=eye(3);
    for i=idx:ref_idx-1
        T=T*i_To_iPlusOne_Transform{i};
    end
elseif idx>ref_idx
    T=eye(3);
    for i=idx-1:ref_idx
        T=T*pinv(i_To_iPlusOne_Transform{i});
    end
else
    T=eye(3);
end

end
