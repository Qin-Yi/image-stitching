% scale space feature detector based upon difference of gaussian filters.
% selects features based upon their maximum response in scale space
function [features,keys] = refine_features(img, pyr, scl, imp,pts, radius)
    
[~,wo]=size(img);
[~,w2]=size(imp{2});

ndx = 1;

% loop through each level of pyramid
for j=2:length(imp)-1                                               

    p=pts{j};                                                      
    img=imp{j};                                                   

    [dh,dw]=size(img);                  
    for i=1:size(p,1)        
        ptx = p(i,2);
        pty = p(i,1);
    
        % ensure neighborhood is not outside of image
        if p(i,1) < radius+3                                      
            p(i,1) = radius+3;
        end
        if p(i,1) > dh-radius-3
            p(i,1) = dh-radius-3;
        end
        if p(i,2) < radius+3
            p(i,2) = radius+3;
        end
        if p(i,2) > dw-radius-3
            p(i,2) = dw-radius-3;
        end
            
        % adjust to sub pixel maxima location
        r = pyr{j}(pty-1:pty+1,ptx-1:ptx+1);                     
        % find center of paraboloid fit to points
        [pcy,pcx] = fit_paraboloid(r);                             
        
        % ignore extreme offsets due to singularities in parabola fitting
        if abs(pcy)>1 || abs(pcx)>1                                              
            pcy=0; pcx=0;
        end
               
        p(i,1) = p(i,1) + pcy;         
        p(i,2) = p(i,2) + pcx;
        ptx = p(i,2);
        pty = p(i,1);
        
        % calculate point locations at pyramid level 2
        px=(pts{j}(i,2)+pcx - 1)*scl^(j-2) + 1;                     
        py=(pts{j}(i,1)+pcy - 1)*scl^(j-2) + 1;
                 
        % calculate Sub-Scale level adjustment
        y1 = interp(pyr{j-1},(p(i,2)-1)*scl+1, (p(i,1)-1)*scl+1);  
        y3 = interp(pyr{j+1},(p(i,2)-1)/scl+1, (p(i,1)-1)/scl+1);
        y2 = interp(pyr{j},p(i,2),p(i,1));
        
        % fit neighborhood of 3 scale points to parabola 
        coef = fit_parabola(0,1,2,y1,y2,y3);  
        % find max in scale space 
        scale_ctr = -coef(2)/2/coef(1);                            
        
        % ignore extreme values due to singularities in parabola fitting
        if abs(scale_ctr-1)>1                                       
            scale_ctr=0;
        end

        % eliminate edge points and enforce minimum separation
        % adjust radius size to account for new scale value
        rad2 = radius * scl^(scale_ctr-1);                          
        % create ring of points at radius around test point
        o=0:pi/8:2*pi-pi/8;                                         
        rx = (rad2)*cos(o);
        ry = (rad2)*sin(o);
        
        rmax = 1e-9;                                                
        rmin = 1e9;
        % get response at feature center
        pval = interp(pyr{j},ptx,pty);                              
        rtst = [];
                
        % check points on ring around feature point
        for k=1:length(rx)   
            % get ring point value with bilinear interpolation
            rtst(k) = interp(pyr{j},ptx+rx(k),pty+ry(k));          
            % calculate distance from feature point for each point in ring
            if pval> 0                                             
                rtst(k) = pval - rtst(k);
            else
                rtst(k) = rtst(k) - pval;
            end
            % test for valid maxima above noise floor                                                       
            rmax = max(rmax,rtst(k));
            rmin = min(rmin,rtst(k));
        end        
        
        % calculate size offset due to edge effects of downsampling
        fac = scl/(wo*2/size(pyr{2},2));                                                               

        % save x and y position IN ORIGINAL IMAGE (sub pixel adjusted)
        features(ndx,1) = (px-1)*wo/w2*fac +1;                      
        features(ndx,2) = (py-1)*wo/w2*fac +1;
        % save scale value (sub scale adjusted in units of pyramid level)
        features(ndx,3) = j+scale_ctr-1;     
        % save x and y position in the current pyramid level
        % i.e. the center pixel is imp{round(features(ndx,3))}(py,px)
        features(ndx,4) = ptx;                                      
        features(ndx,5) = pty;           

        keys(ndx,:) = construct_key(ptx,pty,imp{j},3 * scl^(scale_ctr-1));
        ndx = ndx + 1;
    end
    
end


function v = interp(img,xc,yc)                                      
       % bilinear interpolation between points
       px = round(xc);
       py = round(yc);
       alpha = xc-px;
       beta = yc-py;

       nw = img(py,px);
       ne = img(py,px+1);
       sw = img(py+1,px);
       se = img(py+1,px+1);
       
       v = (1-alpha)*(1-beta)*nw + (1-alpha)*beta*sw + ...
              alpha*(1-beta)*ne + alpha*beta*se;
        
          
 function key = construct_key(px, py, img, sz)
    pct = .75;  
    [h,w] = size(img);    
    [y_off,x_off] = meshgrid(-1:1,-1:1);
    y_off = y_off(:)*pct;
    x_off = x_off(:)*pct;
    
    for i = 1:size(y_off,1)
        % method using interpolated values
        ctrx = px + x_off(i)*sz*2;  
        ctry = py + y_off(i)*sz*2;
        [y,x] = meshgrid(ctry-sz:sz/3:ctry+sz,ctrx-sz:sz/3:ctrx+sz);
        y=y(:); x=x(:);
        t = 0; c = 0;
        for k=1:size(y,1)
            if x(k)<w-1 && x(k)>1 && y(k)<h-1 && y(k)>1 
                t = t + interp(img,x(k),y(k));
                c=c+1;
            end
        end
        t = t/c;      
        key(i) = t;
    end
    key = key/sum(key);
 
    
function [yctr,xctr] = fit_paraboloid(data)

    %Create Solution Matrix
    [yp,xp] = meshgrid(-1:1,-1:1);
    yp=yp(:); xp=xp(:);

    M=zeros(9,6);
    for i=1:length(yp)
        x=xp(i); y=yp(i);
        M(i,:) = [x^2 y^2 x*y x y 1];
    end

    g = data';
    g = g(:);

    cf = (M'*M)\(M'*g);
    a=cf(1); b=cf(2); c=cf(3); d=cf(4); e=cf(5);
    xctr = (-2*b*d+e*c)/(4*a*b-c*c);
    yctr = (-2*a*e-d*c)/(4*a*b-c*c);
    
function [c]=fit_parabola(x1, x2, x3, y1, y2, y3)
    z = [x1^2 x1 1; x2^2 x2 1; x3^2 x3 1];
    c = z\[y1; y2; y3];
    

          



    