function descriptors = SIFTDescriptor(imp, keypoint, Scales)

% SIFTDescriptor Build SIFT descriptors from image at detected key points'
% location with detected key points' scale and angle


    %% Precompute the gradients at all pixels of all pyramid scales
    grad_theta = cell(length(imp),1);
    grad_mag = cell(length(imp),1);
    
    for scale = 1:length(imp)
        currentImage = imp{scale};

        % gradient image
        img_dx = filter2([-1 0 1],currentImage);
        img_dy = filter2([-1;0;1],currentImage);

        % Calculate the magnitude and orientation of the gradient.
        grad_mag{scale} = sqrt(img_dx.^2+img_dy.^2);
        grad_theta{scale} = atan2(img_dy,img_dx);
        % atan2 gives angles from -pi to pi, change that to 0 to 2*pi.
        grad_theta{scale} = mod(grad_theta{scale}, 2*pi);
        
    end
    
    num_angles = 8;
    num_histograms = 4;
    pixelsPerHistogram = 4;
    patch_size = num_histograms * pixelsPerHistogram;
    N = size(keypoint, 1);
    descriptors = zeros(N, num_histograms * num_histograms * num_angles);
            
    for i = 1 : N
        scale = round(Scales(i));    
        % Find the window of pixels that contributes to the descriptor for the current keypoint.
        % center of the DoG keypoint
        xAtScale = keypoint(i, 1);
        yAtScale = keypoint(i, 2);
        x_low = round(xAtScale - patch_size / 2);
        x_high = x_low+patch_size-1;
        y_low = round(yAtScale - patch_size / 2);
        y_high = y_low+patch_size-1;                
            
        % These are the gradient magnitude and angle images from the correct scale level.
        magnitudes = grad_mag{scale};
        thetas = grad_theta{scale};
        try    
            % Extract the patch from that window around the keypoint
            patch_mag = magnitudes(y_low:y_high,x_low:x_high);
            patch_theta = thetas(y_low:y_high,x_low:x_high);
        catch err
            % If any keypoint is too close to the boundary of the image, just skip it.
            continue;
        end
                                                                                     
        % Express gradient directions relative to the dominant gradient direction     
        patch_angle_offset = ComputeDominantDirection(patch_mag, patch_theta);
        patch_theta = patch_theta - patch_angle_offset;
        % re-map patch_theta into the range 0 to 2*pi
        patch_theta = mod(patch_theta, 2*pi);
        % Weight the gradient magnitudes using a gaussian function
        patch_mag = patch_mag .* fspecial('gaussian', patch_size, patch_size / 2);
               
        % Compute the gradient histograms and concatenate them in the          
        % feature variable to form a size 1x128 SIFT descriptor for this keypoint.    
        feature = [];
        for y=1:4:13
            for x=1:4:13
                % For simplicity,just assign all gradient pixels within a square to the same histogram,
                % instead of smoothing a gradient across nearby histograms like SIFT paper
                subdivided_patch_theta = patch_theta(y:y+3,x:x+3);
                subdivided_patch_mag = patch_mag(y:y+3,x:x+3);
                histogram = ComputeGradientHistogram(num_angles,subdivided_patch_mag,subdivided_patch_theta);
                feature=[feature,histogram];
            end
        end
        descriptors(i, :) = feature;
    end
    descriptors = NormalizeDescriptors(descriptors);
end



function histogram = ComputeGradientHistogram(num_bins, gradient_magnitudes, gradient_angles)
% Compute a gradient histogram using gradient magnitudes and directions.
% Each point is assigned to one of num_bins depending on its gradient
% direction; the gradient magnitude of that point is added to its bin.
    angle_step = 2 * pi / num_bins;
    histogram = zeros(1, num_bins);
    [rows,cols] = size(gradient_magnitudes);
    for m=1:rows
        for n=1:cols
            for angle=0:angle_step:(2*pi-angle_step)
                if(angle<=gradient_angles(m,n)&&gradient_angles(m,n)<angle+angle_step)
                    histogram(round(angle/angle_step)+1) = ...
                        histogram(round(angle/angle_step)+1)+gradient_magnitudes(m,n);
                    break;
                end
            end
        end
    end
end


%% Computes the dominant gradient direction for the region around a keypoint
% given the scale of the keypoint and the gradient magnitudes and gradient
% angles of the pixels in the region surrounding the keypoint.
function direction = ComputeDominantDirection(gradient_magnitudes, gradient_angles)
    num_bins = 36;
    histogram = ComputeGradientHistogram(num_bins, gradient_magnitudes, gradient_angles);
    peak=max(histogram);
    loc=find(histogram==peak);
    direction =2*pi*(loc(1)-1)/num_bins;
end


%% normalize all descriptors so they become unit vectors
function descriptors = NormalizeDescriptors(descriptors)
    lengths = sqrt(sum(descriptors.^2, 2));
    lengths(lengths == 0) = 1;
    descriptors = descriptors ./ repmat(lengths, [1 size(descriptors,2)]);

    % suppress large entries
    descriptors(descriptors > 0.2) = 0.2;

    % finally, renormalize to unit length
    lengths = sqrt(sum(descriptors.^2, 2));
    lengths(lengths == 0) = 1;
    descriptors = descriptors ./ repmat(lengths, [1 size(descriptors,2)]);

end