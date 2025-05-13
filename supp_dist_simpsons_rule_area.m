

%{
    
    APPROXIMATES AN AREA UNDER THE CURVE USING SIMPSON'S RULE

%}


function A = supp_dist_simpsons_rule_area(y, h)

    m = length(y);
    
    switch m
        case 1
            A = h * y(1);
        case 2
            A = (h/2) * (y(1) + y(2));
        case 3
            A = (h/3) * (y(1) + 4 * y(2) + y(3));
        case 4
            A = (h/2) * (y(1) + y(2)) + (h/3) * (y(2) + 4 * y(3) + y(4));
        otherwise
            if mod(m,2) == 1
                n = (m-1)/2;
                A = sr(y, h, n);
            else
                n = (m/2) - 1;
                A = (h/2) * (y(1) + y(2)) + sr(y(2:m), h, n);
            end
    end
end

            
function A = sr(y, h, n)            
   A = (h/3) * ( y(1) + 4 * sum(y(2:2:(2*n))) + 2 * sum(y(3:2:(2*n-1))) + y(2*n+1));
end


% -------------------------------------
% Test function...
% -------------------------------------

function try_sr()
    for m = 1:20
        n = m + 1;
        h = (pi/2) / m;
        for i = 1:n
            y(i) = sin((i-1) * pi/(2*m));
        end
        A = fov_simpsons_rule_area(y,h);
        disp(['n = ' num2str(n) ': Area = ' num2str(A, '%f')]);
    end
end
