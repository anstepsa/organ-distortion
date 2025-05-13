



function supp_dist_visualise_volume(NEW_FIG, SMOOTH, v, d)

    if NEW_FIG
        hf = figure; set(hf, 'Color','w');
    end
    
    v = double(v);

    if SMOOTH, v = smooth3(v,'gaussian',15,2); end

    v_dim = size(v);
    
    [x,y,z] = meshgrid(1:v_dim(2),1:v_dim(1),1:v_dim(3));
    x = x * d(1); y = y * d(2); z = z * d(3);
  
    is = isosurface(y,x,z,v,0.01);
    p = patch(is);
    
    set(p,'FaceColor','red','EdgeColor','none', 'BackFaceLighting','reverselit');
    xlabel('y'); ylabel('x'); zlabel('z');
    view(3); axis xy
    camlight; lighting phong
    rotate3d; daspect([1 1 1]);
    
    if NEW_FIG
        close(hf);
    end
    
end
