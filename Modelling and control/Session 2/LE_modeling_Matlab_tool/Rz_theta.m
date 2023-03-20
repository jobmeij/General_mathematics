% function Rz = Rz_theta(th, Matr_dimension)

function Rz = Rz_theta(th, Matr_dimension)

if Matr_dimension == 3,
    Rz = [cos(th) -sin(th) 0
          sin(th)  cos(th) 0
                0        0 1];
    
elseif Matr_dimension == 2,
    Rz = [cos(th) -sin(th)
          sin(th)  cos(th)];
end