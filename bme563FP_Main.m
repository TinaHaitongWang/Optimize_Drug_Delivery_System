% BME 563 Final Project 
% Author: Haitong Wang 

% Main file 
clear all; close all; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable Amax, V_L, t_hat 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Amax = [45, 50, 55]; % cm2, half space 
V_L = [1,1.5,2]; % mL = cm3 
t_hat = [1,1.5,2]; % hrs
V = [1.5,4]; % gel volume 

% v_opt = zeros(7,3,3);
% 
% % find the optimal intial volume that achieve A(t_hat) = Amax, maximize the
% % average SF values over all the variations 
% v0 = [1:0.01:2]; 
% for i = 1:length(Amax)
%     indexVolume = i
%     for j = 1: length(t_hat)
%        indexTime = j
%        tSpan = 0:0.25:(t_hat(j)*3600); % 2 hr interval 
%        v_opt(:,i,j) = findVolumeOpt(Amax(i),tSpan,v0); 
%     end 
% end  
% 
% save('v_opt.mat','v_opt')

load('v_opt.mat') % loading the optimal volume calculated in the previous section
% gel order: 3000, 3001, 3002, 4002, DG1, DG2, DG3

% adjust the initial gel volume
V = [1.5,4.00]; % gel volume 
for i = 1:length(v_opt(:,1,1))
    for j = 1: length(v_opt(1,:,1))
        for k = 1:length(v_opt(1,1,:))
            if(v_opt(i,j,k)<V(1))
                v_opt(i,j,k) = V(1);
            elseif (v_opt(i,j,k)>V(2))
                v_opt(i,j,k) = V(2);
            end
        end
    end
end

% calculate A(t) and h(t) for each combination using the optimal volume
% three time span
Amax = [45, 50, 55]; % cm2, half space
V_L = [1,1.5,2]; % mL = cm3
t_hat = [1,1.5,2]; % hrs
t_1 = 0:0.25:(t_hat(1)*3600);
t_2 = 0:0.25:(t_hat(2)*3600);
t_3 = 0:0.25:(t_hat(3)*3600);
h0 = 0.1 ; %cm
w = 2; % cm
D = 6E-6 ; % cm^2/s
F_4ml = 4.44822; %  1 lbf = 4.44822 N =
k = F_4ml/4;
F = k.*v_opt;

% Gel 3000
m = 630*0.1*(1/10^2)^2; % 0.1 Pa*s^n = N*s^n/ m^2 = N*s^n/ cm^2
tau_0 = 2*(1/10^2)^2; % 1 Pa = N/m^2 = N/m^2 *(1 m^2/)
n = .455;
for i = 1:length(Amax)
    [A_3000,h_3000] = calculateA_yieldstress(m,tau_0,n,t_1,v_opt(1,i,1),F(1,i,1));
    A_3000t1(:,i)= A_3000;  h_3000t1(:,i) = h_3000;
    [A_3000,h_3000] = calculateA_yieldstress(m,tau_0,n,t_2,v_opt(1,i,2),F(1,i,2));
    A_3000t2(:,i)= A_3000;  h_3000t2(:,i) = h_3000;
    [A_3000,h_3000] = calculateA_yieldstress(m,tau_0,n,t_3,v_opt(1,i,3),F(1,i,3));
    A_3000t3(:,i)= A_3000;  h_3000t3(:,i) = h_3000;
end

% Gel 3001
m = 254*0.1*(1/10^2)^2; % 0.1 Pa*s^n = N*s^n/ m^2 = N*s^n/ cm^2
n = .569;
for i = 1:length(Amax)
    [A_3001,h_3001 ] = calculateA_withoutYS(m,n,t_1,v_opt(2,i,1),F(2,i,1));
    A_3001t1(:,i)= A_3001;  h_3001t1(:,i) = h_3001;
    [A_3001,h_3001 ]  = calculateA_withoutYS(m,n,t_2,v_opt(2,i,2),F(2,i,2));
    A_3001t2(:,i)= A_3001;  h_3001t2(:,i) = h_3001;
    [A_3001,h_3001 ]  = calculateA_withoutYS(m,n,t_3,v_opt(2,i,3),F(2,i,3));
    A_3001t3(:,i)= A_3001;  h_3001t3(:,i) = h_3001;
end

% Gel 3002
m = 484*0.1*(1/10^2)^2; % 0.1 Pa*s^n = N*s^n/ m^2 = N*s^n/ cm^2
n = .518;
for i = 1:length(Amax)
    [A_3002,h_3002] = calculateA_withoutYS(m,n,t_1,v_opt(3,i,1),F(3,i,1));
    A_3002t1(:,i)= A_3002;  h_3002t1(:,i) = h_3002;
    [A_3002,h_3002] = calculateA_withoutYS(m,n,t_2,v_opt(3,i,2),F(3,i,2));
    A_3002t2(:,i)= A_3002;  h_3002t2(:,i) = h_3002;
    [A_3002,h_3002] = calculateA_withoutYS(m,n,t_3,v_opt(3,i,3),F(3,i,3));
    A_3002t3(:,i)= A_3002;  h_3002t3(:,i) = h_3002;
end

% Gel 4002
m = 816*0.1*(1/10^2)^2; % 0.1 Pa*s^n = N*s^n/ m^2 = N*s^n/ cm^2
tau_0 = 20*(1/10^2)^2; % 1 Pa = N/m^2 = N/m^2 *(1 m^2/)
n = .309;
for i = 1:length(Amax)
    [A_4002,h_4002] = calculateA_yieldstress(m,tau_0,n,t_1,v_opt(4,i,1),F(4,i,1));
    A_4002t1(:,i)= A_4002;  h_4002t1(:,i) = h_4002;
    [A_4002,h_4002] = calculateA_yieldstress(m,tau_0,n,t_2,v_opt(4,i,2),F(4,i,2));
    A_4002t2(:,i)= A_4002;  h_4002t2(:,i) = h_4002;
    [A_4002,h_4002] = calculateA_yieldstress(m,tau_0,n,t_3,v_opt(4,i,3),F(4,i,3));
    A_4002t3(:,i)= A_4002;  h_4002t3(:,i) = h_4002;
end
% Gel DG1
m = 662*0.1*(1/10^2)^2; % 0.1 Pa*s^n = N*s^n/ m^2 = N*s^n/ cm^2
tau_0 = 2*(1/10^2)^2; % 1 Pa = N/m^2 = N/m^2 *(1 m^2/)
n = .512;
for i = 1:length(Amax)
    [A_DG1,h_DG1] = calculateA_yieldstress(m,tau_0,n,t_1,v_opt(5,i,1),F(5,i,1));
    A_DG1t1(:,i)= A_DG1;  h_DG1t1(:,i) = h_DG1;
    [A_DG1,h_DG1] = calculateA_yieldstress(m,tau_0,n,t_2,v_opt(5,i,2),F(5,i,2));
    A_DG1t2(:,i)= A_DG1;  h_DG1t2(:,i) = h_DG1;
    [A_DG1,h_DG1] = calculateA_yieldstress(m,tau_0,n,t_3,v_opt(5,i,3),F(5,i,3));
    A_DG1t3(:,i)= A_DG1;  h_DG1t3(:,i) = h_DG1;
end

% Gel DG2
m = 928*0.1*(1/10^2)^2; % 0.1 Pa*s^n = N*s^n/ m^2 = N*s^n/ cm^2
tau_0 = 38*(1/10^2)^2; % 1 Pa = N/m^2 = N/m^2 *(1 m^2/)
n = .450;
for i = 1:length(Amax)
    [A_DG2,h_DG2] = calculateA_yieldstress(m,tau_0,n,t_1,v_opt(6,i,1),F(6,i,1));
    A_DG2t1(:,i)= A_DG2;  h_DG2t1(:,i) = h_DG2;
    [A_DG2,h_DG2] = calculateA_yieldstress(m,tau_0,n,t_2,v_opt(6,i,2),F(6,i,2));
    A_DG2t2(:,i)= A_DG2;  h_DG2t2(:,i) = h_DG2;
    [A_DG2,h_DG2] = calculateA_yieldstress(m,tau_0,n,t_3,v_opt(6,i,3),F(6,i,3));
    A_DG2t3(:,i)= A_DG2;  h_DG2t3(:,i) = h_DG2;
end
% Gel DG3
m = 57*0.1*(1/10^2)^2; % 0.1 Pa*s^n = N*s^n/ m^2 = N*s^n/ cm^2
n = .618;
for i = 1:length(Amax)
    [A_DG3,h_DG3] = calculateA_withoutYS(m,n,t_1,v_opt(7,i,1),F(7,i,1));
    A_DG3t1(:,i)= A_DG3;  h_DG3t1(:,i) = h_DG3;
    [A_DG3,h_DG3] = calculateA_withoutYS(m,n,t_2,v_opt(7,i,2),F(7,i,2));
    A_DG3t2(:,i)= A_DG3;  h_DG3t2(:,i) = h_DG3;
    [A_DG3,h_DG3] = calculateA_withoutYS(m,n,t_3,v_opt(7,i,3),F(7,i,3));
    A_DG3t3(:,i)= A_DG3;  h_DG3t3(:,i) = h_DG3;
end
 
% % compute the M(t) for each combination  
for i = 1: length(Amax)
    M3000t1(i,:) = ComputeMt(A_3000t1(:,i),h_3000t1(:,i),t_1,Amax(i),v_opt(1,i,1));
    M3000t2(i,:) = ComputeMt(A_3000t2(:,i),h_3000t2(:,i),t_2,Amax(i),v_opt(1,i,2));
    M3000t3(i,:) = ComputeMt(A_3000t3(:,i),h_3000t3(:,i),t_3,Amax(i),v_opt(1,i,3));
    M3001t1(i,:) = ComputeMt(A_3001t1(:,i),h_3001t1(:,i),t_1,Amax(i),v_opt(2,i,1));
    M3001t2(i,:) = ComputeMt(A_3001t2(:,i),h_3001t2(:,i),t_2,Amax(i),v_opt(2,i,2));
    M3001t3(i,:) = ComputeMt(A_3001t3(:,i),h_3001t3(:,i),t_3,Amax(i),v_opt(2,i,3));
    M3002t1(i,:) = ComputeMt(A_3002t1(:,i),h_3002t1(:,i),t_1,Amax(i),v_opt(3,i,1));
    M3002t2(i,:) = ComputeMt(A_3002t2(:,i),h_3002t2(:,i),t_2,Amax(i),v_opt(3,i,2));
    M3002t3(i,:) = ComputeMt(A_3002t3(:,i),h_3002t3(:,i),t_3,Amax(i),v_opt(3,i,3));
    M4002t1(i,:) = ComputeMt(A_4002t1(:,i),h_4002t1(:,i),t_1,Amax(i),v_opt(4,i,1));
    M4002t2(i,:) = ComputeMt(A_4002t2(:,i),h_4002t2(:,i),t_2,Amax(i),v_opt(4,i,2));
    M4002t3(i,:) = ComputeMt(A_4002t3(:,i),h_4002t3(:,i),t_3,Amax(i),v_opt(4,i,3));
    MDG1t1(i,:) = ComputeMt(A_DG1t1(:,i),h_DG1t1(:,i),t_1,Amax(i),v_opt(5,i,1));
    MDG1t2(i,:) = ComputeMt(A_DG1t2(:,i),h_DG1t2(:,i),t_2,Amax(i),v_opt(5,i,2));
    MDG1t3(i,:) = ComputeMt(A_DG1t3(:,i),h_DG1t3(:,i),t_3,Amax(i),v_opt(5,i,3));
    MDG2t1(i,:) = ComputeMt(A_DG2t1(:,i),h_DG2t1(:,i),t_1,Amax(i),v_opt(6,i,1));
    MDG2t2(i,:) = ComputeMt(A_DG2t2(:,i),h_DG2t2(:,i),t_2,Amax(i),v_opt(6,i,2));
    MDG2t3(i,:) = ComputeMt(A_DG2t3(:,i),h_DG2t3(:,i),t_3,Amax(i),v_opt(6,i,3));
    MDG3t1(i,:) = ComputeMt(A_DG3t1(:,i),h_DG3t1(:,i),t_1,Amax(i),v_opt(7,i,1));
    MDG3t2(i,:) = ComputeMt(A_DG3t2(:,i),h_DG3t2(:,i),t_2,Amax(i),v_opt(7,i,2));
    MDG3t3(i,:) = ComputeMt(A_DG3t3(:,i),h_DG3t3(:,i),t_3,Amax(i),v_opt(7,i,3));
end

% compute the G(t) for each combination  
for i = 1: length(Amax)
    for j = 1: length(V_L)
        Gt3000t1(i,j) = ComputeGt(A_3000t1(end,i),Amax(i),V_L(j),v_opt(1,i,1));
        Gt3000t2(i,j) = ComputeGt(A_3000t2(end,i),Amax(i),V_L(j),v_opt(1,i,2));
        Gt3000t3(i,j) = ComputeGt(A_3000t3(end,i),Amax(i),V_L(j),v_opt(1,i,3));
        Gt3001t1(i,j) = ComputeGt(A_3001t1(end,i),Amax(i),V_L(j),v_opt(2,i,1));
        Gt3001t2(i,j) = ComputeGt(A_3001t2(end,i),Amax(i),V_L(j),v_opt(2,i,2));
        Gt3001t3(i,j) = ComputeGt(A_3001t3(end,i),Amax(i),V_L(j),v_opt(2,i,3));
        Gt3002t1(i,j) = ComputeGt(A_3002t1(end,i),Amax(i),V_L(j),v_opt(3,i,1));
        Gt3002t2(i,j) = ComputeGt(A_3002t2(end,i),Amax(i),V_L(j),v_opt(3,i,2));
        Gt3002t3(i,j) = ComputeGt(A_3002t3(end,i),Amax(i),V_L(j),v_opt(3,i,3));
        Gt4002t1(i,j) = ComputeGt(A_4002t1(end,i),Amax(i),V_L(j),v_opt(4,i,1));
        Gt4002t2(i,j) = ComputeGt(A_4002t2(end,i),Amax(i),V_L(j),v_opt(4,i,2));
        Gt4002t3(i,j) = ComputeGt(A_4002t3(end,i),Amax(i),V_L(j),v_opt(4,i,3));
        GtDG1t1(i,j) = ComputeGt(A_DG1t1(end,i),Amax(i),V_L(j),v_opt(5,i,1));
        GtDG1t2(i,j) = ComputeGt(A_DG1t2(end,i),Amax(i),V_L(j),v_opt(5,i,2));
        GtDG1t3(i,j) = ComputeGt(A_DG1t3(end,i),Amax(i),V_L(j),v_opt(5,i,3));
        GtDG2t1(i,j,:) = ComputeGt(A_DG2t1(end,i),Amax(i),V_L(j),v_opt(6,i,1));
        GtDG2t2(i,j,:) = ComputeGt(A_DG2t2(end,i),Amax(i),V_L(j),v_opt(6,i,2));
        GtDG2t3(i,j,:) = ComputeGt(A_DG2t3(end,i),Amax(i),V_L(j),v_opt(6,i,3));
        GtDG3t1(i,j,:) = ComputeGt(A_DG3t1(end,i),Amax(i),V_L(j),v_opt(7,i,1));
        GtDG3t2(i,j,:) = ComputeGt(A_DG3t2(end,i),Amax(i),V_L(j),v_opt(7,i,2));
        GtDG3t3(i,j,:) = ComputeGt(A_DG3t3(end,i),Amax(i),V_L(j),v_opt(7,i,3));
    end 
end 

% compute the SF(t) for 27 combination 
SF_3000 = zeros(3,3,3); 
SF_3001 = zeros(3,3,3);
SF_3002 = zeros(3,3,3);
SF_4002 = zeros(3,3,3);
SF_DG1 = zeros(3,3,3);
SF_DG2 = zeros(3,3,3);
SF_DG3 = zeros(3,3,3);

for i = 1: 3 % Amax
    for j = 1: 3 % VL 
        SF_3000(i,j,1) = M3000t1(i).*Gt3000t1(i,j);
        SF_3000(i,j,2) = M3000t2(i).*Gt3000t2(i,j);
        SF_3000(i,j,3) = M3000t3(i).*Gt3000t3(i,j);
        SF_3001(i,j,1) = M3001t1(i).*Gt3001t1(i,j);
        SF_3001(i,j,2) = M3001t2(i).*Gt3001t2(i,j);
        SF_3001(i,j,3) = M3001t3(i).*Gt3001t3(i,j);
        SF_3002(i,j,1) = M3002t1(i).*Gt3002t1(i,j);
        SF_3002(i,j,2) = M3002t2(i).*Gt3002t2(i,j);
        SF_3002(i,j,3) = M3002t3(i).*Gt3002t3(i,j);
        SF_4002(i,j,1) = M4002t1(i).*Gt4002t1(i,j);
        SF_4002(i,j,2) = M4002t2(i).*Gt4002t2(i,j);
        SF_4002(i,j,3) = M4002t3(i).*Gt4002t3(i,j);
        SF_DG1(i,j,1) = MDG1t1(i).*GtDG1t1(i,j);
        SF_DG1(i,j,2) = MDG1t2(i).*GtDG1t2(i,j);
        SF_DG1(i,j,3) = MDG1t3(i).*GtDG1t3(i,j);
        SF_DG2(i,j,1) = MDG2t1(i).*GtDG2t1(i,j);
        SF_DG2(i,j,2) = MDG2t2(i).*GtDG2t2(i,j);
        SF_DG2(i,j,3) = MDG2t3(i).*GtDG2t3(i,j);
        SF_DG3(i,j,1) = MDG3t1(i).*GtDG3t1(i,j);
        SF_DG3(i,j,2) = MDG3t2(i).*GtDG3t2(i,j);
        SF_DG3(i,j,3) = MDG3t3(i).*GtDG3t3(i,j);
    end 
end 

% find the optimal SF value and v_opt value for each gel 
SF_hat = zeros(7,1);
v_hat = zeros(7,1); 
SF_error = zeros(7,3);
SF_value = zeros(27,7);
% gel order: 3000, 3001, 3002, 4002, DG1, DG2, DG3
% if different initial volume produce the same SF value, we want to pich
% the one has the smallest leackage volume, longer time, and smaller
% initial volume 
[SF_hat(1), v_hat(1),SF_error(1,:),SF_value(:,1)] = ComputeOptimalValue(SF_3000,v_opt,1);
[SF_hat(2), v_hat(2),SF_error(2,:),SF_value(:,2)] = ComputeOptimalValue(SF_3001,v_opt,2);
[SF_hat(3), v_hat(3),SF_error(3,:),SF_value(:,3)] = ComputeOptimalValue(SF_3002,v_opt,3);
[SF_hat(4), v_hat(4),SF_error(4,:),SF_value(:,4)] = ComputeOptimalValue(SF_4002,v_opt,4);
[SF_hat(5), v_hat(5),SF_error(5,:),SF_value(:,5)] = ComputeOptimalValue(SF_DG1,v_opt,5);
[SF_hat(6), v_hat(6),SF_error(6,:),SF_value(:,6)] = ComputeOptimalValue(SF_DG2,v_opt,6);
[SF_hat(7), v_hat(7),SF_error(7,:),SF_value(:,7)] = ComputeOptimalValue(SF_DG3,v_opt,7);


