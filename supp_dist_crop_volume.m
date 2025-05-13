
%{

    Crop a volume and associated spatial reference object...

%}
function [Vc, Rc, ok] = supp_dist_crop_volume(V, R, margin)
    
    Vc = [];
    Rc = [];
    
    stats = regionprops(V, 'BoundingBox');
    BB = stats.BoundingBox;

    ul = BB(1:2); ul = ul - margin; 
    
    if any(ul < 1),       ok = -1; disp('Margin too big!'); return; end
    if ul(1) > size(V,1), ok = -1; disp('Margin too big!'); return; end
    if ul(2) > size(V,2), ok = -1; disp('Margin too big!'); return; end
    
    wi = BB(4:5); wi = wi + 2 * margin;
    
    if any((ul + wi) < 1),          ok = -1; disp('Margin too big!'); return; end
    if (ul(1) + wi(1)) > size(V,1), ok = -1; disp('Margin too big!'); return; end
    if (ul(2) + wi(2)) > size(V,2), ok = -1; disp('Margin too big!'); return; end
    
    xmin = ul(1); xmax = ul(1) + wi(1);
    ymin = ul(2); ymax = ul(2) + wi(2);

    limits = [xmin, xmax, ymin, ymax, NaN, NaN]; 

    Vc = subvolume(V, limits);

    zWorldLimits = R.ZWorldLimits;
    
    xWorldLimits = [R.XWorldLimits(1) + (xmin-1) * R.PixelExtentInWorldX, ...
                    R.XWorldLimits(1) + (xmax-1) * R.PixelExtentInWorldX];
    
    yWorldLimits = [R.YWorldLimits(1) + (ymin-1) * R.PixelExtentInWorldY, ...
                    R.YWorldLimits(1) + (ymax-1) * R.PixelExtentInWorldY];

    Rc = imref3d(size(Vc), xWorldLimits, yWorldLimits, zWorldLimits);
    
    disp('Volume cropped!');  % <abg> 22-05-2024
    
    ok = 1;

end


