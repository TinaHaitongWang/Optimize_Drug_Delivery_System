% BME 563 Final Project
% Author: Haitong Wang

% function that compute volume for each combination of Amax and t_hat

function [v_opt] = findVolumeOpt(Amax,t_hat,v0)
% gel order: 3000, 3001, 3002, 4002, DG1, DG2, DG3
v_opt = zeros(1,7);
for i = 1:7
    currentGel = i
    v_index = 0;
    v = v0;
    k = 0;
    v_new =0;
    while(true)
        k = k +1; 
        if (k ==1)
        [A,H]= estimateA(v0,t_hat,i);
        end 
        for j = 1: length(A)
            if (abs(A(j)-Amax)<=0.01) ||  (abs(Amax-A(j))<=0.01)
                v_index = j;
                v_opt(i) = v(v_index);
            end
        end
        
        if (v_index == 0)
            if (A(1)>=Amax)            
                v_new = v(1);
                v = [(v_new-0.5):0.001:(v_new+0.1)];
                [A,H] = estimateA(v,t_hat,i);
            elseif (A(end)<=Amax)
                v_new = v(end);
                v = [(v_new-0.1):0.001:(v_new+0.5)];
                [A,H] = estimateA(v,t_hat,i);
            else
                for j = 1: length(A)
                    if (fix(A(j)) == Amax)
                        if (v_new == v(j))
                            v_new = v(j);
                            v = [(v_new-0.001):0.00005:(v_new+0.001)];
                        else 
                        v_new = v(j);
                        v = [(v_new-0.01):0.0005:(v_new+0.01)];
                        end 
                        [A,H] = estimateA(v,t_hat,i);
                        break;
                    end
                end
            end
        else
            break;
        end
        
    end
    
end

end