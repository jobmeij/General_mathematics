% function Ry = Ry_theta(th, Matr_dimension)

function Ry = Ry_theta(th, Matr_dimension)

if Matr_dimension == 3,
    Ry = [cos(th) 0 sin(th)
                0 1       0
         -sin(th) 0 cos(th)];
    
elseif Matr_dimension == 2,
    Ry = [ cos(th)  sin(th)
          -sin(th)  cos(th)];
end