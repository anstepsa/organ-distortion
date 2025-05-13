




% In-built function 'rms' is in Statistics Toolbox...
function r = my_rms(a)
    try
        r = rms(a);
    catch
        r = sqrt(mean(a.^2));
    end
end