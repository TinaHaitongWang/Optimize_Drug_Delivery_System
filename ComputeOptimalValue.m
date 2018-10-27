% BME 563 Final Project 
% Author: Haitong Wang 

% function to find the optimal SF value and correponding v value 
function [SF_hat, v_hat,SF_error,SF_value] = ComputeOptimalValue(SF,v_opt,gel)
SF_value = reshape(SF,[27,1]);
pd = fitdist(SF_value,'Normal');
SF_hat = pd.mu;
ci = paramci(pd);
error = pd.mu - ci(1,1);
SF_error = [ci(1,1),ci(2,1),error];
v = v_opt(gel,:,:);
vv = v(:);
v_hat = mean(vv);
end 