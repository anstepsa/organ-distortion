


function I_t = my_translate_by(I, T, R)
% Translate an image in real world by 'T' (mm) keeping in same spatial reference object...

    if nargin < 1
        R.ImageSize = [128,128];
        R.PixelExtentInWorldX = 1;
        R.PixelExtentInWorldY = 1;
        T = [20,40];
        imshow(zeros(128,128));
        h = impoly();
        pos = getPosition(h);
        I = poly2mask(pos(:,1), pos(:,2), 128, 128);
    end

    figure
    imshow(I,[]);


    R2d = imref2d(R.ImageSize(1:2), R.PixelExtentInWorldX, R.PixelExtentInWorldY);

    T2d  = affine2d([ 1,    0,  0;
                      0,    1,  0;
                    T(1), T(2), 1]);
               
    I_t = imwarp(I, R2d, T2d, 'OutputView',R2d);           

    figure
    imshow(I_t);
    
end