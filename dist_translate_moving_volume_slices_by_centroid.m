


function [tV1, tV1C, translate_by] = dist_translate_moving_volume_slices_by_centroid(V1, V2, R, SURFACE_THRESHOLD)
 
    % Find centroid of each slice of V1...
    V1C = ones(2,size(V1,3)) * -1;
    for slice = 1:size(V1,3)
        stats = regionprops((V1(:,:,slice) > SURFACE_THRESHOLD),'centroid');
        if ~isempty(stats)
            V1C(:,slice) = stats.Centroid; % Centroid is returned by regionprops as (X,Y)
        end
    end
    
    % Find centroid of each slice of V2...
    V2C = ones(2,size(V2,3)) * -1;
    for slice = 1:size(V2,3)
        stats = regionprops((V2(:,:,slice) > SURFACE_THRESHOLD),'centroid');
        if ~isempty(stats)
            V2C(:,slice) = stats.Centroid; 
        end
    end

    % Slice spatial reference object...
    % Is 'sR1', were it to be formed, identical...?
    sR2 = imref2d([size(V2,1),size(V2,2)], R.XWorldLimits, R.YWorldLimits);

    % Initialise variables...
    translate_by = NaN(2,size(V1,3));
    tV1          = zeros(size(V2));
    tV1C         = ones(2,size(V1,3)) * -1;
    
    for slice = 1:size(V1,3)
        
        if all(V1(:,:,slice) < SURFACE_THRESHOLD)
            continue
        end
        if all(V1C(:,slice) == -1) || all(V2C(:,slice) == -1)
            continue
        end
        
        % Determine translation to take centroid to centroid in slice)...
        translate_by(:,slice) = (V2C(:,slice) - V1C(:,slice)) .* [R.PixelExtentInWorldX, R.PixelExtentInWorldY]';
        
        % Turn this into a transformation...
        trans = affine2d([1,                     0,                     0;
                          0,                     1,                     0;
                          translate_by(1,slice), translate_by(2,slice), 1]);
                      
        % Translate the moving slice by this amount - so that the slice centroids
        % should now coincide...
        tV1(:,:,slice) = imwarp(V1(:,:,slice), sR2, trans, 'OutputView',sR2); % WHY NOT 'sR1' or equivalent here (first use)...?

        % Get the new centroids for moving slices...
        stats = regionprops((tV1(:,:,slice) > SURFACE_THRESHOLD),'centroid');
        if ~isempty(stats)
            tV1C(:,slice) = stats.Centroid; 
        end
        
    end
    
end
