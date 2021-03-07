function [area_vec] = CalcFaceAreas(Solver_setup,idx)

%---------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------
%-------------------Calculate Area of faces of geometry to determine how large------------------------
%---------------------------------------------------------------------------------------------
triangle_areas = Solver_setup.triangle_area_m2;
    for n = 1:(length(idx)-1)
        area_vec(n) = sum(triangle_areas(idx(n):idx(n+1)))
    end   
    

end

