function x_min=AllocateDataFun(V, unitCost, slope_m, tau, fog_unprocessed, fog_processed, fog_capacity, arrival_rate, DC_que, unitBandCost, bandWidth)
sec2hr = 3600;%unit: from second to hour
w2mw = 10^6; %unit: from W to MW
w = 0.1;
beta  = 0.01;
num_DC = size(DC_que,2);

R=10000;

% the first parameter represents the amount of workload assigned to edge device
coeffient = zeros(num_DC,num_DC*2+1);
for j=1:num_DC
    coeffient(j,j+1) = 1;
    coeffient(j,j+1+num_DC) = 1;
end
% A = [ones(1,num_DC+1) zeros(1,num_DC) 
%     zeros(1,num_DC+1) ones(1,num_DC) 
%     coeffient];
A = [R ones(1,num_DC) zeros(1,num_DC) 
    zeros(1,num_DC+1) ones(1,num_DC) 
    coeffient];
b = [fog_unprocessed fog_processed bandWidth];
lb = zeros(1,num_DC*2+1);
ub = [min(fog_capacity*tau,arrival_rate), bandWidth, bandWidth];

factor_1 = V*unitCost*slope_m/w2mw/sec2hr + w*fog_processed - fog_unprocessed;
factor_2 = V*unitBandCost + DC_que - fog_unprocessed;
factor_3 = V*unitBandCost + beta*DC_que - fog_processed;
f = [factor_1,factor_2,factor_3];
options = optimoptions('linprog','Algorithm','dual-simplex','Display','off'); 
%options = optimoptions('linprog','Algorithm','interior-point','Display','off'); 
%options = optimoptions('linprog','Algorithm','dual-simplex'); 
x_min = linprog(f,A,b,[],[],lb,ub,options);
end