

function [V1_reg, V2_interp, R_full, t_transl] = dist_register_volumes(v1, R1, v2, R2, PL)
% Do a rigid body translation to register 'moving' volume v1 in space R1 to 'fixed' volume v2 in space R2...
% Express the result in space 'R_full' encompassing both R1 and R2...
% Interpolate fixed volume 'v2' to this new space too...

    [optimizer, metric] = imregconfig('monomodal');

    % Find the translation which accomplishes the registration...
    t_transl = imregtform(v1, R1, v2, R2, 'translation', optimizer, metric, 'PyramidLevels', PL, 'DisplayOptimization', true);
    
    % Note that this translation can be spurious (in absolute terms) if the
    % ImagePositionPatient DICOM header field is wrong in some way (e.g. 
    % registered or re-formatted images). This has no effect on subsequent 
    % distortion calculations i.e. once the translation has been applied.
    
    % Construct a spatial reference object (R_full) which encompasses both R1 and R2...
    New.XWorldLimits = [min(R1.XWorldLimits(1), R2.XWorldLimits(1)), max(R1.XWorldLimits(2), R2.XWorldLimits(2))];
    New.YWorldLimits = [min(R1.YWorldLimits(1), R2.YWorldLimits(1)), max(R1.YWorldLimits(2), R2.YWorldLimits(2))];
    New.ZWorldLimits = [min(R1.ZWorldLimits(1), R2.ZWorldLimits(1)), max(R1.ZWorldLimits(2), R2.ZWorldLimits(2))];
    
    New.X_dim = round((New.XWorldLimits(2) - New.XWorldLimits(1)) / R2.PixelExtentInWorldX);
    New.Y_dim = round((New.YWorldLimits(2) - New.YWorldLimits(1)) / R2.PixelExtentInWorldY);
    New.Z_dim = round((New.ZWorldLimits(2) - New.ZWorldLimits(1)) / R2.PixelExtentInWorldZ);
    
    R_full = imref3d([New.Y_dim, New.X_dim, New.Z_dim], New.XWorldLimits, New.YWorldLimits, New.ZWorldLimits);
    
    
    % If the combined spatial reference object is unreasonably big (this can
    % happen with oblique acquisitions re-formatted to axial where the 
    % ImagePosition becomes ?corrupted), just use the fixed volume space
    %(for speed and accuracy of subsequent calculations)...
    
    % [Could possibly introduce error if moving volume outside reference
    % space but in context this is unlikely.]
    
    if any( (R_full.ImageSize ./ R2.ImageSize) > 1.5)
        warning('Combined spatial reference is large');%: reverting to fixed volume reference space.');
       % R_full = R2;
    end
     if any( (R_full.ImageSize(3) ./ R2.ImageSize(3)) > 1.5)
        warning('Large spatial reference offset in ''z''!');%: reverting to fixed volume reference space.');
       % R_full = R2;
     end
    
    % Translate in 3-D the moving volume giving result in 'full' space...
    [V1_reg,    R_full] = imwarp(v1, R1, t_transl,         'OutputView',R_full);
    
    % Interpolate fixed 'v2' volume into this space too (i.e. warp by identity transform)...
    [V2_interp, R_full] = imwarp(v2, R2, affine3d(eye(4)), 'OutputView',R_full);
    
end

