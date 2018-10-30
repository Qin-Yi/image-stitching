% scale space feature detector based upon difference of gaussian filters.
% selects features based upon their maximum response in scale space
function local = find_features(pyr,scl,threshold,radius)
levels = size(pyr);
levels = levels(2);
    
for i=2:levels-1
    % find maxima and minima at current scale level
    mx = find_peak(pyr{i},threshold,radius); 
    mn = find_peak(-pyr{i},threshold,radius);
    
    % find coords in level above 
    mx2 = round((mx-1)/scl) + 1;  
    mn2 = round((mn-1)/scl) + 1;  
    % do neighbor comparison in scale space above
    mx_above = comp_neighbor(pyr{i},pyr{i+1},mx,mx2,radius);   
    mn_above = comp_neighbor(-pyr{i},-pyr{i+1},mn,mn2,radius);

    if i>1
        % find coords in level below
        mx2 = round((mx-1)*scl) + 1; 
        mn2 = round((mn-1)*scl) + 1;
        % do comparison in scale below
        mx_below = comp_neighbor(pyr{i},pyr{i-1},mx,mx2,radius);  
        mn_below = comp_neighbor(-pyr{i},-pyr{i-1},mn,mn2,radius);
        
        % get coord list for retained maxima and minima
        maximum = mx(find(mx_below & mx_above),:); 
        minimum = mn(find(mn_below & mn_above),:);    
    else
        maximum = mx(find(mx_above),:);
        minimum = mn(find(mn_above),:);
    end
    
    % combine maxima and minima into list for return
    local{i} = [maximum; minimum];                                   
end

%   Compare a vector of pixels with its neighbors in another scale 
function v = comp_neighbor(img1,img2,i,i2,radius)                    
    % i and i2 are column vectors of r,c coords
    
    if (size(i2,1))==0 || size(img2,1)<11 || size(img2,2)<11
        v=zeros(length(i),1);
    else
        
        [h,w] = size(img1);
        [h2,w2] = size(img2);
    
        % create set of offsets within radius 
        [y,x]=meshgrid(-20:20,-20:20);                              
        z = (x.^2+y.^2)<=radius^2;
        [y,x]=find(z);
        x=x-21; y=y-21;
    
        % create boundary listing
        bound = ones(size(i2,1),2)*[h2-radius 0;0 w2-radius];        
        % test bounds to make all points within image
        i2 = i2 - ((i2 > bound).*(i2-bound+1));                      
        i2 = i2 + ((i2 < radius+1).*(radius-i2+1));
        
        % create indices from x,y coords
        i2 = vec(i2,h2);                                             
        i = vec(i,h);
    
        p = img1(i);
        res = ones(length(i),1);
    
        % check against all points within radius
        for j=1:length(x)                                            
            itest = i2 + x(j) + h2*y(j);
            p2 = img2(itest);
            res = res & (p>=p2);
        end

        % store results in binary vector
        v = res;                                                     
    end
    
function v = vec(points,h)
    y = points(:,1);
    x = points(:,2);    
    v = y + (x-1)*h;  %create index vectors
    
    
%% finds local maxima within a grayscale image.  
% Each point is checked against all of the pixels within a given radius to be a local max/min.
% The magnitude of pixel values must be above the given threshold to be picked as a valid maxima or minima. 
function [mx] = find_peak(img,threshold,radius)
[h,w] = size(img);

% get interior image subtracting radius pixels from border
p = img(radius+1:h-radius, radius+1:w-radius);  

% get pixels above threshold
[yp,xp] = find(p>threshold);     
yp = yp+radius; xp = xp+radius;
pts = yp+(xp-1)*h;

% create offset list for immediate neighborhood
z=ones(3,3);                    
z(2,2)=0;
[y,x]=find(z);
y=y-2; x=x-2;

if size(pts,2)>size(pts,1)
    pts = pts';
end

% test for max within immediate neighborhood
if size(pts,1)>0
    maxima=ones(length(pts),1);
    for i=1:length(x)
        pts2 = pts + y(i) + x(i)*h;
        maxima = maxima & img(pts)>img(pts2);
    end
    % create new index list of good points
    xp = xp(find(maxima)); 
    yp = yp(find(maxima));
    pts = yp+(xp-1)*h;                              
end
    
% create offset list for radius of pixels
[y,x]=meshgrid(-20:20,-20:20);  
% include points within radius without immediate neighborhood
z = (x.^2+y.^2)<=radius^2 & (x.^2+y.^2)>(radius)^2;   
[y,x]=find(z);
x=x-21; y=y-21;

maxima = ones(length(pts),1);

% test within radius of pixels (done after first test for slight speed increase)
if size(pts,1)>0
    for i = 1:length(x)
        pts2 = pts + y(i) + x(i)*h;
        maxima = maxima & img(pts)>img(pts2);   
    end
    xp = xp(find(maxima));                
    yp = yp(find(maxima));                            
    mx = [yp xp];    
else
    mx = [];
end
    

    