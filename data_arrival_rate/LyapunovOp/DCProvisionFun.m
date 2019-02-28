function [num, speed]=DCProvisionFun(DC_que, V, unitCost, rho, power_min, speed_option, density, theta, num_min, num_max)
sec2hr = 3600;%unit: from second to hour
minVal = +inf;
num_speed_option = size(speed_option,2);
bit2GB = 8*2^30; %unit: from bit to MB
w2mw = 10^6; %unit: from W to MW
ghz2hz = 10^9;  %unit: from GHz to Hz
for j=1:num_speed_option
    inc_speed = (speed_option(j) - speed_option(1));
    f = V*unitCost*(rho*inc_speed^theta + power_min)/w2mw/sec2hr - DC_que*speed_option(j)*ghz2hz/(bit2GB*density);
    intcon = 1;
    lb = num_min;
    ub = num_max;
    %no display the detail information of the optimization process
    options = optimoptions('intlinprog','Display','off');
    [x,fval] = intlinprog(f,intcon,[],[],[],[],lb,ub,options);
    if(fval<minVal)
        minVal = fval;
        minNum = x;
        minSpeed = speed_option(j);
    end
end
num = minNum;
speed = minSpeed;
end