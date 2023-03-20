% function Rx = Rx_theta(th, Matr_dimension)

function Rx = Rx_theta(th, Matr_dimension)

if Matr_dimension == 3,
    Rx = [1       0        0
          0 cos(th) -sin(th)
          0 sin(th)  cos(th)];
    
elseif Matr_dimension == 2,
    Rx = [ cos(th)  -sin(th)
           sin(th)  cos(th)];
end