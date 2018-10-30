function H = RANSACFit(p1, p2, dscpt1, dscpt2)
% Use RANSAC to find a robust transformation 

    match = SIFTSimpleMatcher(dscpt1, dscpt2);
    N = size(match, 1);
    if N<4
        error('not enough matches to produce a transformation matrix')
    end
    num_iter = 200;
    seedSetSize = max(ceil(0.2 * N),3);
    maxInlierError = 30;
    goodFitThresh = floor(0.7 * N);
    H = eye(3);
    
    a = Inf;
    for i = 1 : num_iter
        [sample, rest] = random(match, seedSetSize);
        h = homography(p1(sample(:, 1), :), p2(sample(:, 2), :));
        err = ComputeError(h, p1, p2, rest);
        flag = (err <= maxInlierError);
        if sum(flag(:)) + seedSetSize >= goodFitThresh
            point = [sample; rest(flag, :)];
            h = homography(p1(point(:, 1), :), p2(point(:, 2), :));
            b = sum(ComputeError(h, p1, p2, point));
            if b < a
                H = h;
                a = b;
            end
        end
    end
end

function dists = ComputeError(H, pt1, pt2, match)
% Compute the error using transformation matrix H to transform the point in pt1 to its matching point in pt2.
% Error is measured as the Euclidean distance between (transformed pt1) and pt2 in homogeneous coordinates.
    transform_pt1 = H*[pt1(match(:,1),:)';ones(1,size(match,1))];
    subtract = pt2(match(:,2),:)-transform_pt1(1:2,:)';
    dists = sqrt(subtract(:,1).^2+subtract(:,2).^2);
end

function [D1, D2] = random(D, splitSize)
    idx = randperm(size(D, 1));
    D1 = D(idx(1:splitSize), :);
    D2 = D(idx(splitSize+1:end), :);
end

function match = SIFTSimpleMatcher(descriptor1, descriptor2)
% Each descriptor from descriptor1 can at most be matched to one member of descriptor2, 
% but descriptors from descriptor2 can be matched more than once.

threshold = 0.7;
match = [];

[N1,~] = size(descriptor1);
[N2,~] = size(descriptor2);
for i=1:N1
    %   For each descriptor vector in descriptor1, find the Euclidean distance
    %   between it and each descriptor vector in descriptor2.
    distance=[];
    for j=1:N2
        subtract = descriptor1(i,:)-descriptor2(j,:);
        distance = [distance,norm(subtract)];
    end
    % If the smallest distance is less than thresh*(the next smallest distance), 
    % we say that the two vectors are a match
    sort_distance = sort(distance);
    if (sort_distance(1) < threshold*sort_distance(2))
        j = find(distance==sort_distance(1));
        match=[match;[i,j]];
    end
end

end

function H = homography( Pt1, Pt2 )
%   Computes the transformation matrix that transforms a point from coordinate frame 1 to coordinate frame 2
    N = size(Pt1,1);
    if size(Pt1, 1) ~= size(Pt2, 1)
        error('Dimensions unmatched.');
    elseif N<4
        error('At least 4 points are required.');
    end
    
    x_ref = Pt1(:,1);
    y_ref = Pt1(:,2);
    x_src = Pt2(:,1);
    y_src = Pt2(:,2);
    
    A = zeros(N*2,8);
    
    A(1:2:end,1:3) = [x_ref, y_ref, ones(N,1)];
    A(2:2:end,4:6) = [x_ref, y_ref, ones(N,1)];
    A(1:2:end,7:8) = [-x_ref.*x_src, -y_ref.*x_src];
    A(2:2:end,7:8) = [-x_ref.*y_src, -y_ref.*y_src];

    B = [x_src, y_src];
    B = reshape(B',N*2,1);

    h = A\B;
    H = [h(1),h(2),h(3);h(4),h(5),h(6);h(7),h(8),1];
end







