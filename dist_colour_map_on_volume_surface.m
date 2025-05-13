


function hf = dist_colour_map_on_volume_surface(config, params, d_array, Ts, V, R)
% Show the distortion mapped by colour onto the volume surface...
% 
%   d_array : an array of local distortion values indexed by [phi, slice]
%   Ts      : an array of slice translations mapping centroid of moving volume slice to centroid of fixed volume slice
%   V       : the volume on which to show the surface colour map    
%   R       : spatial reference object associated with volume above
%
    if nargin < 1
        % Optional line to load variables...
        load('temp.mat', 'config', 'params', 'd_array', 'Ts', 'V1_reg', 'R_full', 'results');
        V = V1_reg; R = R_full; 
    end
 
    % Get parameters...
    SURFACE_THRESHOLD     = config.SURFACE_THRESHOLD;  
    
    SMOOTH_IT             = params.SMOOTH_VOLUME_FOR_DISPLAY;
    DISTORTION_MAP_MAX_MM = params.DISTORTION_MAP_MAX_MM;
    
    % Get the slice range without zero-filled slices at top and bottom...
    valid_slices = get_non_zero_filled_slices(V, SURFACE_THRESHOLD);
    
    % Crop the data to these slices only...
    V_crop = V(:, :, valid_slices);
    d_array_crop = d_array(:, valid_slices);
    
    % Get rid of translation data for non-valid slices...
    Ts_crop = Ts(:,valid_slices);
    
    % Find the new centroid of this volume...
    stats = regionprops((V_crop > SURFACE_THRESHOLD),'centroid');
    VC_crop = stats.Centroid; 
    
    % Smooth the final image (for cosmetic effect only)
    if SMOOTH_IT
        V_crop = supp_dist_imgaussian(V_crop, 5.0);
    end
    
    % Get the distortion field as a Cartesian matrix from the distortion
    % array (polar coordinates by slice)...
    distortion_field = get_distortion_field(V_crop, VC_crop, d_array_crop, Ts_crop, R); % [r,c,s]

    % N.B. Volume is in coordinates system [rows, cols, Z], so we must
    % produce a grid in [rows, cols, Z] as well...
    
    
    % Meshgrid (in 'mm')... N.B. 'ndgrid' output has different dimensions to meshgrid
    %[X,Y,Z] = ndgrid( (1:size(V_crop,1)) * R.PixelExtentInWorldY, (1:size(V_crop,2)) * R.PixelExtentInWorldX, (1:size(V_crop,3)) * R.PixelExtentInWorldZ);
    [rows_Y, cols_X, Z] = ndgrid( (1:size(V_crop,1)) * R.PixelExtentInWorldY, (1:size(V_crop,2)) * R.PixelExtentInWorldX, (1:size(V_crop,3)) * R.PixelExtentInWorldZ);
    
    % Calculate the surface to plot...in space with coordinates (rows, cols, Z)...
    
    %[faces, verts, colors] = isosurface(X,Y,Z, V_crop, 0.5, distortion_field);
    [faces, verts, colors] = isosurface(rows_Y, cols_X, Z, V_crop, 0.5, distortion_field);

    % Show the figure...
    hf = figure; set(hf, 'Color','w');
    
    % Construct a colour map...
    colormap(gca, cat(1, flip(jet(100),1), [0.7,0.7,0.7]));  % Pre-append 'grey' to the jet map
    
    % Scale colour map...
    set(gca, 'CLim',[-DISTORTION_MAP_MAX_MM, DISTORTION_MAP_MAX_MM]); % Colormap indexed from [-max,max] 'mm distortion'
    
    % Set figure parameters...
    set(gca, 'DataAspectRatio', [1, 1, 1]);

    % Apply lighting...
    camlight;
    
    % Show the surface...
    patch('Vertices',verts, ...
          'Faces',faces, ... 
          'FaceVertexCData',colors, ... 
          'FaceColor','interp', ...           
          'EdgeColor','interp', ...
          'CDataMapping','scaled', ...
          'FaceLighting','phong', ...
          'EdgeLighting','phong', ...
          'AmbientStrength',0.8);% ...
    
    % Prepare the 'tree' of centroid displacement vectors...  
    translate_slice(1,:) = Ts_crop(1,~isnan(Ts_crop(1,:)));
    translate_slice(2,:) = Ts_crop(2,~isnan(Ts_crop(2,:)));
    
    % ... and plot it...
    hold on;
    ht = line([0,0], [0,0], [0,size(translate_slice,2)] .* R.PixelExtentInWorldZ); 
    set(ht, 'LineStyle','-', 'LineWidth',1.0, 'Color','k');
    for s = 1:size(translate_slice,2)
        
        % To be drawn in space [rows, cols, Z]... translate_slice is in
        % [cols_Xs, rows_Ys] so must be flipped...
        %my_arrow(gca, [0,translate_slice(1,s)], [0,translate_slice(2,s)], [s,s] .* R.PixelExtentInWorldZ);   
        my_arrow(gca, [0,translate_slice(2,s)], [0,translate_slice(1,s)], [s,s] .* R.PixelExtentInWorldZ);           
        
    end
    hold off;
    
    % Final adjustments to figure...
    axis vis3d; 
    
    rotate3d;
    view([-15,20]);
    
    ylabel('R -> L (mm)');
    xlabel('P -> A (mm)');
    zlabel('I -> S (mm)');
    
    colorbar;
    
    drawnow;
    
end


function my_arrow(ax, Xs, Ys, Zs)
% Shows a displacement vector in axes 'ax'...

    if     (Xs(1) == 0 && isnan(Xs(2))) ...
        && (Ys(1) == 0 && isnan(Ys(2)))
        return
    end
    
    ht = line(Xs, Ys, Zs); 
    set(ht, 'LineStyle','-', 'LineWidth',2.0, 'Color','r');
    hp = plot3(ax, Xs(2), Ys(2), Zs(2));
    set(hp, 'Marker','o', 'MarkerSize',5, 'MarkerFaceColor','r', 'MarkerEdgeColor','none');
    
end


function D = get_distortion_field(V, VC, d_array, T, R)
% Returns a 'distortion field' (D) associated with the volume V in space r and distortion
% array 'd_array'...

    % Fill 'not calculated' values as +999 - mapped to end of colormap i.e.
    % will appear grey in final plot...
    d_array(isnan(d_array)) = 999;

    ni = size(V,2); % i,  X,  cols
    nj = size(V,1); % j,  Y,  rows
    
    slices = size(V,3);

    [X, Y] = meshgrid(1:ni, 1:nj);  % X is in form [rows, cols]
                                    % Y is in form [rows, cols]
    X = X - VC(1);
    Y = Y - VC(2);
    
    D = zeros(nj, ni, slices); % [rows, cols, slice]
    
    ps = 0;
    
    for s = 1:slices
        for i = 1:ni
            for j = 1:nj

                % Convert polar array 'd_array' to values at volume points
                % in Cartesian array...
                
                phi = round(atand(abs( Y(j,i) / X(j,i) ))); 

                if X(j,i) == 0 && Y(j,i) >= 0, phi =  90; end
                if X(j,i) == 0 && Y(j,i) <= 0, phi = 270; end
                if Y(j,i) == 0 && X(j,i) >= 0, phi =   0; end
                if Y(j,i) == 0 && X(j,i) <= 0, phi = 180; end

                if phi == 0, phi = 1; end

                if X(j,i) > 0 && Y(j,i) < 0, phi = 360 - phi; end
                if X(j,i) < 0 && Y(j,i) < 0, phi = 180 + phi; end
                if X(j,i) < 0 && Y(j,i) > 0, phi = 180 - phi; end

                if isnan(phi), phi = 1; end

                h = s;
                
                if s > ps
                    ps = s;
                    disp(['Processing slice: ' num2str(ps)]);
                end
                
                try
                    % Construct the distortion field element...
                    D(j, i, s) = d_array(phi, h);   % (rows, cols, slice)
                catch ME
                    disp(ME.message);
                end

            end
        end
        
        % Apply slice centroid translation to local distortion field (this
        % gives correct colouring to surface - which is on the 'original
        % object')...
        
        if ~any(isnan(T(:,s)))
            
            % D is [r,c,Z]; T is [cols_Xs, rows_Ys] as distance between
            % centroids...
            D(:,:,s) = my_translate_by(D(:,:,s), -T(:,s), R); 
            
        end
        
    end

end



function I_t = my_translate_by(I, T, R)
% Translate an image in real world by [Tx_cols, Ty_rows] (mm) keeping in same spatial reference object...

    R2d = imref2d(R.ImageSize(1:2), R.PixelExtentInWorldX, R.PixelExtentInWorldY);

    T2d  = affine2d([ 1,    0,  0;
                      0,    1,  0;
                    T(1), T(2), 1]);
               
    I_t = imwarp(I, R2d, T2d, 'OutputView',R2d);           

end




function slices = get_non_zero_filled_slices(V, SURFACE_THRESHOLD)

    start = -1; last = -1;

    for s = 1:size(V,3)
        if ~isempty(nonzeros(V(:,:,s) > SURFACE_THRESHOLD))
            if start < 0, start = s; break; end
        end
    end
    for s = size(V,3):-1:1
        if ~isempty(nonzeros(V(:,:,s) > SURFACE_THRESHOLD))
            if last < 0, last = s; break; end
        end
    end

    slices = start:last;
    
end
