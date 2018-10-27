% BME 563 Final Project
% Author: Haitong Wang

% function that compute Area for estimated initial volume at each t_hat 

function [A,h] = estimateA(v,t_hat,gel)

    F_4ml = 4.44822; 
    k = F_4ml/(4 );
    F = k.*v;
    A = zeros(length(v),1);
    h = zeros(length(v),1);
    
    for i = 1:length(v)
        if (gel ==1)
            % Gel 3000
            m = 630*0.1*(1/10^2)^2; % 0.1 Pa*s^n = N*s^n/ m^2 = N*s^n/ cm^2
            tau_0 = 2*(1/10^2)^2; % 1 Pa = N/m^2 = N/m^2 *(1 m^2/)
            n = .455;
            [A_3000,h_3000 ] = calculateA_yieldstress(m,tau_0,n,t_hat,v(i),F(i));
            A(i) = A_3000(end);
            h(i) = h_3000(end);
        elseif (gel ==2)
            % Gel 3001
            m = 254*0.1*(1/10^2)^2; % 0.1 Pa*s^n = N*s^n/ m^2 = N*s^n/ cm^2
            n = .569;
            [A_3001,h_3001 ] = calculateA_withoutYS(m,n,t_hat,v(i),F(i));
            A(i) = A_3001(end);
            h(i) = h_3001(end);
        elseif (gel ==3)
            % Gel 3002
            m = 484*0.1*(1/10^2)^2; % 0.1 Pa*s^n = N*s^n/ m^2 = N*s^n/ cm^2
            n = .518;
            [A_3002,h_3002 ] = calculateA_withoutYS(m,n,t_hat,v(i),F(i));
            A(i) = A_3002(end);
            h(i) = h_3002(end);
        elseif (gel ==4)
            % Gel 4002
            m = 816*0.1*(1/10^2)^2; % 0.1 Pa*s^n = N*s^n/ m^2 = N*s^n/ cm^2
            tau_0 = 20*(1/10^2)^2; % 1 Pa = N/m^2 = N/m^2 *(1 m^2/)
            n = .309;
            [A_4002,h_4002] = calculateA_yieldstress(m,tau_0,n,t_hat,v(i),F(i));
            A(i) = A_4002(end);
            h(i) = h_4002(end);
        elseif (gel ==5)
            % Gel DG1
            m = 662*0.1*(1/10^2)^2; % 0.1 Pa*s^n = N*s^n/ m^2 = N*s^n/ cm^2
            tau_0 = 2*(1/10^2)^2; % 1 Pa = N/m^2 = N/m^2 *(1 m^2/)
            n = .512;
            [A_DG1,h_DG1] = calculateA_yieldstress(m,tau_0,n,t_hat,v(i),F(i));
            A(i) = A_DG1(end);
            h(i) = h_DG1(end);
        elseif (gel ==6)
            % Gel DG2
            m = 928*0.1*(1/10^2)^2; % 0.1 Pa*s^n = N*s^n/ m^2 = N*s^n/ cm^2
            tau_0 = 38*(1/10^2)^2; % 1 Pa = N/m^2 = N/m^2 *(1 m^2/)
            n = .450;
            [A_DG2,h_DG2] = calculateA_yieldstress(m,tau_0,n,t_hat,v(i),F(i));
            A(i) = A_DG2(end);
            h(i) = h_DG2(end);
        elseif (gel ==7)
            % Gel DG3
            m = 57*0.1*(1/10^2)^2; % 0.1 Pa*s^n = N*s^n/ m^2 = N*s^n/ cm^2
            n = .618;
            [A_DG3,h_DG3] = calculateA_withoutYS(m,n,t_hat,v(i),F(i));
            A(i) = A_DG3(end);
            h(i) = h_DG3(end);
        end
        
    end

end

