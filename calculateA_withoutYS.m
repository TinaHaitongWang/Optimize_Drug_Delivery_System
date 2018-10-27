% BME 563 Final Project 
% Author: Haitong Wang 

% function calculating spreading A for gel without yield stress 
function [A,h_t] = calculateA_withoutYS(m,n,t_f,v,F)
w = 2; %cm 
h0 = 0.1 ; %cm


p = (2*n+3)/n;
G = (((2*n+1)/(4*n))^n) * (m/(8*n+16)) * ((v^(n+2))/(w^(n+1)));

h_t = (((1/h0).^p)+ (p*((F/G).^(1/n)).*t_f)).^(-1/p);
h_t = h_t';
A = v ./ h_t;

end 