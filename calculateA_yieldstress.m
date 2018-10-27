% BME 563 Final Project 
% Author: Haitong Wang 

% function calculating spreading A for gel with yield stress 
% A is half space 

function  [A,h_t] = calculateA_yieldstress(m,tau_0,n,tSpan,v,F)
w = 2; %cm 

C = (2*w*m)/(n+2)*((2*n+1)/n)^n*(v/(4*w))^(n+2);
Q = tau_0*(v^2/(16*w));
 
y0 = 0.1 ; %cm 
[t,y] = ode45(@(t,y) odefcn(t,y,n,C,Q,F),tSpan, y0);

A = v ./ y;
h_t = y;

    function dydt = odefcn(t,y,n,C,Q,F)   
          dydt = - C^(-1/n)*y^((3*n+2)/n)*(F*y-Q*(y^-2))^(1/n);
    end 

end 