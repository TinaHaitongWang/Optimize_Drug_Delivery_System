% BME 563 Final Project 
% Author: Haitong Wang 

% Compute the M(t) function 
function M = ComputeMt(A,h,tspan,Amax,v)
D = 6E-6 ; % cm^2/s

z = 1./(h(:).^2); 
soln_z(:) = cumtrapz(tspan,z);

for i = 2: (length(soln_z(:))) % each time point
    e = exp(-((pi/2)^2) * D * soln_z(i));
    for j = 1 : 700 % n from 0 to infinite
        e = e + exp(-(((j*pi)+(pi/2))^2)*D*soln_z(i));
    end
    e_t(i) = e;
end

soln_Mt = zeros(1,length(A(:)));
for i = 1:length(A(:))% time
    if(A(i)<Amax)
        soln_Mt(i) = A(i)/h(i);
    elseif (A(i)>=Amax)
        soln_Mt(i) = Amax/h(i);
    end
end
y = 2*D.*soln_Mt(2:end).* e_t(2:end);
Mt(:) = cumtrapz(tspan(2:end),y);
M= Mt(end)/v;

end 