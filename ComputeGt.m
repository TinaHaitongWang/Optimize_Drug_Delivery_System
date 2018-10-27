% BME 563 Final Project 
% Author: Haitong Wang 

% Compute the G(t) function

function Gt = ComputeGt(A,Amax,V_L,v)
if (V_L>=v)
    A_star = 1000000;
else
    A_star = Amax/(1-(V_L/v));
end
    
if(A<=Amax)
    Gt = A/Amax;
elseif (A>=Amax) && (A<=A_star)
    Gt = (A_star - A)/(A_star-Amax);
elseif( A_star < A)
    Gt = 0;
end

end 