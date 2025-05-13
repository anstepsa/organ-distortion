




function [d, sem_d, d_res, sem_d_res, d_array, d_array_res] = dist_compute_radial_distortion_by_slice(config, V1, V1C, V2, R, Ts)

    % ----------------------------------------------------------
    % Initialise...
    phi_range = 1:360;
    phi_range = phi_range * pi / 180.0;  % ...in radians
    
    d_rho     = config.RADIAL_STEP;
    d_r       = config.RADIAL_STEP;
    
    n_rho     = (1 / d_rho) * config.SEARCH_RADIUS_PIX;
    
    rho_range = 0:d_rho:(n_rho * d_rho);

    np        = length(phi_range) * size(V1,3);
    
    d_array            = NaN(length(phi_range), size(V1,3));
    d_array_res        = NaN(length(phi_range), size(V1,3));
    
    dist_magnitude     = NaN(np,1); 
    dist_magnitude_res = NaN(np,1); 
    
    s = 1;
    % ----------------------------------------------------------
    
    % For each slice in moving volume...
    for slice = 1:size(V1,3)
        
        % If centroid of moving volume slice NOT in fixed volume...
        if ~(is_in_region(V2(:,:,slice), V1C(:,slice)'))
            warning('Translated centroid of moving slice not in reference field of view...');
            continue;  
        end
        
        % For all radial projections...
        for phi = phi_range
            
            at_boundary  = false;
            abort_search = false;
            
            % Search for boundary in 'moving' volume (V1) starting from origin (slice centroid) along radial
            % projection (distance along this is labelled 'rho')...
            for rho = rho_range
                
                [x, y] = pol2cart(phi, rho);
                
                % Find pixel value of pixel in moving volume (V1) nearest (phi, rho) where
                % origin is at centroid of the slice of V1...
                pix_val = nearest_neighbour(V1(:,:,slice), V1C(:,slice)', [x, y]);
                
                % Are we at a surface boundary of V1...?
                if pix_val < config.SURFACE_THRESHOLD   % I.e. we've found the edge of V1 slice...at 'rho'!
                  
                    if rho == 0
                        warning('rho == 0 : boundary found at origin...');
                        abort_search = true;  
                        %continue   % THIS NEEDS LOOKING AT
                    end

                    % Initialise search displacement radius at value stretching back to origin...
                    r = -rho;
                    
                    % Search along a radial projection for the boundary of
                    % fixed volume (V2): distance along this projection, starting at V1 boundary, is
                    % 'r' and is ranged from '-rho' to 
                    while ~at_boundary && ~abort_search
                        
                        [xd, yd] = pol2cart(phi, rho + r);
                        
                        % Find pixel value of pixel in V2 nearest (rho, phi), where
                        % origin is at centroid of the slice of V1 (not V2)...
                        pix_val = nearest_neighbour(V2(:,:,slice), V1C(:,slice)', [xd, yd]);
                        
                        if pix_val < config.SURFACE_THRESHOLD  % I.e. found boundary of V2...
                            
                            at_boundary = true;
                 
                            % The following are in 'mm'...
                            dist_local         = r  * [ (cos(phi) * R.PixelExtentInWorldX), ...
                                                        (sin(phi) * R.PixelExtentInWorldY) ];

                            dist_resultant     = dist_local + Ts(:,slice)';  

                            dist_magnitude(s)  = sign(r) * sqrt( (r * cos(phi) * R.PixelExtentInWorldX) .^2 + ...
                                                                 (r * sin(phi) * R.PixelExtentInWorldY) .^2 );
                                                                    
                            dist_magnitude_res(s) = sqrt((dist_resultant(1)^2 + dist_resultant(2)^2));

                            d_array(    round(phi * 180.0 / pi), slice) = dist_magnitude(s);
                            d_array_res(round(phi * 180.0 / pi), slice) = dist_magnitude_res(s);

                            if config.VERBOSE > 0                          
                                disp(['s = ' num2str(s) ...
                                      ' --- (phi, rho, slice) = (' num2str(round(phi * 180/pi), '%d') ', ' num2str(rho, '%.2f') ', ' num2str(slice, '%d') ') ' ...
                                      ' --- r = ' num2str(r, '%.2f') ' --- d = ' num2str(dist_magnitude(s), '%.4f') ' --- resultant_d = ' num2str(dist_magnitude_res(s), '%.4f')]);
                            end

                            s = s + 1;
                            
                        end
                        
                        r = r + d_r;
                        
                    end
                    
                end
                
                if at_boundary, break; end
                if abort_search, break; end
                
            end
                        
        end
        
        disp(['-----> [' num2str(slice) ' / ' num2str(size(V1,3)) '] slices completed...']);
        
    end
    
    d     = my_rms(dist_magnitude(    ~isnan(dist_magnitude)    ));
    d_res = my_rms(dist_magnitude_res(~isnan(dist_magnitude_res)));
 
    pix_sz    = R.PixelExtentInWorldX;
    
    sem_d     = pix_sz / sqrt(numel(dist_magnitude(    ~isnan(dist_magnitude)    )));
    sem_d_res = pix_sz / sqrt(numel(dist_magnitude_res(~isnan(dist_magnitude_res))));
 
end


function tf = is_in_region(R, pt)
% Is a point in a slice region 'R'...
% R is an image with non-zero values in the region to be considered...

    val = nearest_neighbour(R, [0,0], pt);
    tf = (val > 0);

end


function val = nearest_neighbour(Vs, V1C, p_in)
% Return the image pixel value at the nearest pixel point in slice image Vs to
% point 'p_in', where p_in has coordinates based on origin 'V1C' and is expressed with
% arbitrary precision...

    % Vs is (row,col)
    % V1C is (x,y)
    % p is (x,y)

    %if V1C(1) == -1 && V1C(2) == -1, val = 0; return; end

    p = p_in + V1C;
    
    p(p < 1) = 1;  % Trap p elements to be at least 1.0... stops bugs in remainder of this function...
    
    v_dim = size(Vs);
    
    if any( (v_dim - [p(2),p(1)] ) <= 0)
        warning('Nearest neighbour might be outside field of view: getting best guess...');
        try
            val = Vs(floor(p(2)),floor(p(1)));
        catch ME
            if floor(p(2)) < 1, val = Vs(1, floor(p(1))); return; end
            if floor(p(1)) < 1, val = Vs(floor(p(2)), 1); return; end
            if floor(p(2)) > v_dim(1), val = Vs(v_dim(1), floor(p(1))); return; end
            if floor(p(1)) > v_dim(2), val = Vs(floor(p(2)), v_dim(2)); return; end
            val = 0;
            disp(ME.message);
        end
        return
    end

    % (The following is much quicker than using 'interp2')
    x_range = floor(p(1)):floor(p(1))+1;
    y_range = floor(p(2)):floor(p(2))+1;
    
    S = Vs(y_range, x_range);

    D = zeros(4,1); count = 1;
    for x = 0:1
        for y = 0:1
            D(count) = (p(1) - x) ^ 2 + (p(2) - y) ^ 2;
            count = count + 1;
        end
    end
    [~,i] = min(D);
    ind = [0,0;  0,1;  1,0;  1,1];
    
    val = S(ind(i,1)+1, ind(i,2)+1);
   
end


