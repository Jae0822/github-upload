%% -----Generate problem instances---------
M = 30;     % number of source nodes and fog gateways
K = 4;     % number of data centers
Iters = 100; %running for one day (1440)
tau = 60;   %unit time of each slot, 100 millisecond.The time interval of data acquisitionn is 60s (1 minute).

Mbps2GBps= 8*2^10; %unit: from bit to MB
power_idle = 325; %unit: W

slope_m = 30.6/8; %unit: Joule per Gigabit to GB/s
fog_capacity = 1.8/8; %unit: Gbps to GB/s

num_min = 20;
num_max = 1000;

speed_option = [2.0 2.5 3.0 3.5 4.0];   %in GHz
ghz2hz = 10^9;  %unit: from GHz to Hz

w2mw = 10^6; %unit: from W to MW

rho = 120;     %the coefficient of power consumption
power_min = 120; %the minimum power consumption required to maintain the server in the active state
theta = 2;      %exponent parameter
sec2hr = 3600;%unit: from second to hour

gamma = 200;  %processing desity of DCs (in cycles/bit)
bit2GB= 8*2^30; %unit: from bit to GB
x_max = tau*num_max*speed_option(end)*ghz2hz/(bit2GB*gamma); %the unit of x_max is GB

w = 0.1;     %the weight of past data used to complex analysis
V = 400;
beta  = 0.01;

%% -----Algorithms-----------------------------------------
% Define the set of data arrival rates
a_options=[1 3 9 14 18 20];

num_a = size(a_options,2);
result_cost = zeros(num_a,1);
QueLen_bar = zeros(num_a,1);
avg_speed = zeros(num_a,1);
avg_num = zeros(num_a,1);

throughput = zeros(num_a,1);
cost_thr = zeros(num_a,1);

sum_bandCost = zeros(num_a,Iters);  % total cost for transmission
sum_migrationCost = zeros(num_a,Iters); %total cost for migration
sum_energy = zeros(num_a,Iters);  %total energy consumption
cost_matrix = zeros(num_a,Iters);
sum_energy_fog = zeros(num_a,Iters);  %total energy consumption

G_avg = zeros(num_a,Iters);
Z_avg = zeros(num_a,Iters);
Q_avg = zeros(num_a,Iters);
QueLen = zeros(num_a,Iters);

%u = 18;
variance = 1;

for a_variable = 1:num_a
    a_index = a_variable;
    a = (variance*randn(Iters,M) + a_options(a_variable));
    s = zeros(Iters,K);     % vector for aggeragator DCs
    mu2 = zeros(Iters,K);       % the amount of data arriving at each DC per slot
    bandCost = zeros(Iters,M);
    fog_bandCost = zeros(Iters,M);
    power = zeros(Iters,K);     %power consumption of each DC
    energy = zeros(Iters,K);     %energy consumption of each DC
    
    num = zeros(Iters,K);   % the number of active servers in each DC
    speed = zeros(Iters,K);
    x = zeros(Iters,K); % processing rate of data center
    r = zeros(Iters,K);     % arrival rate from fog gateways and other DCs
    unitPowerCost = unifrnd(30,70,Iters,K);
    unitPowerCost_fog = unifrnd(30,70,Iters,M);
    
    migDataFromDC = zeros(Iters,K); %the processed data
    migDataFromFog = zeros(Iters,K);
    migCostFromFog = zeros(Iters,K);
    % Original Queues
    G = zeros(Iters,M);
    Z = zeros(Iters,M);
    Q = zeros(Iters,K);
    
    fog_B = zeros(Iters,M);
    power_fog = zeros(Iters,M);
    energy_fog = zeros(Iters,M);
    
    %%-----Initialization-----------------------------------------------
    
    for iter = 1:Iters
        u2=500;
        variance2=50;
        bandWidth = tau*(variance2*randn(M,K) + u2)/Mbps2GBps;
        unitBandCost = 0.05*ones(M,K); % unit cost of different connectivity options (0.05 $/GB)
        
        unitMigrationCost = 0.05*ones(K,K);    %unit cost(0.05 $/GB)
        
        m = zeros(K,K);     % vector for migration data of DCs
        mm = zeros(K,K);     % vector for migration costs of DCs
        
        m_tmp = zeros(K,K); %temp migration data
        mm_tmp = zeros(K,K); %temp migration cost
        preData_tmp = zeros(1,K);
        preCost_tmp = zeros(1,K);
        
        mu = zeros(M,K);    % data allocation matrix
        z = zeros(M,K);    % data allocation matrix
        allocation_matrix = zeros(M,K*2+1);
        
        for i=1:M
            if(a(iter,i)<0)
              a(iter,i) = 0;
            end
            for j=1:K
                if(bandWidth(i,j)<0)
                    bandWidth(i,j)=0;
                end
            end
            
            allocation_matrix(i,:) = AllocateDataFun(V, unitPowerCost_fog(iter,i), slope_m, tau, G(iter,i), Z(iter,i), fog_capacity, a(iter,i), Q(iter,:), unitBandCost(i,:),bandWidth(i,:));
            fog_B(iter,i) = allocation_matrix(i,1);
            mu(i,:) = allocation_matrix(i,2:ceil(end/2));
            z(i,:) = allocation_matrix(i,ceil(end/2)+1:end);
            power_fog(iter,i) = (slope_m*fog_B(iter,i)/tau+power_idle)/w2mw;
            energy_fog(iter,i) = unitPowerCost_fog(iter,i)*power_fog(iter,i)*tau/sec2hr;
            
            bandCost(iter,i) = sum(unitBandCost(i,:).*mu(i,:))+sum(unitBandCost(i,:).*z(i,:));  %select one connectivity option randomly, unit: GB
            
            %----Update Queue------
            G(iter+1,i) = max(G(iter,i) - fog_B(iter,i) - sum(mu(i,:)),0) + a(iter,i);
            Z(iter+1,i) = max(Z(iter,i) - sum(z(i,:)),0) + w*fog_B(iter,i);
        end
        sum_energy_fog(a_index,iter) = sum(energy_fog(iter,:));
        sum_bandCost(a_index,iter) = sum(bandCost(iter,:));
        
        %find the index of aggegrator DC
        for j = 1:K         %assuming DC j is the aggegrator DC
            mu2(iter,j) = sum(mu(:,j))+beta*sum(z(:,j));
            for jj = 1:K
                if((iter>1) && (jj~=j))
                    
                    m_tmp(jj,j) = m_tmp(jj,j) + w*migDataFromDC(iter-1,jj);    %the amount of migration data
                    
                    mm_tmp(jj,j) = unitMigrationCost(jj,j)*m_tmp(jj,j);
                    m_tmp(jj,j) = beta*m_tmp(jj,j);
                else
                    m_tmp(jj,j) = 0;
                    mm_tmp(jj,j) = 0;
                end
            end
        end
        
        [index, s(iter,:)] = SelectDC(Q(iter,:),m_tmp,mm_tmp,V);  % DC selection function
        subscript = index;
        
        y = zeros(K,K);    % data migration vector in different slots
        
        for j = 1:K
            
            [suit_num, suit_speed] = DCProvisionFun(Q(iter,j),V,unitPowerCost(iter,j),rho,power_min,speed_option,gamma,theta,num_min,num_max); %processing rate of of data center
            num(iter,j) = suit_num;
            speed(iter,j) = suit_speed;
            x(iter,j) = tau*num(iter,j)*speed(iter,j)*ghz2hz/(bit2GB*gamma);
            
            %calculate power consumption of DC
            incremental_speed = speed(iter,j) - speed_option(1);
            power(iter,j) = num(iter,j)*(rho*incremental_speed^theta + power_min)/w2mw; % factor*(1-rho) = 120 w; factor = 300 unit: Mw
            
            energy(iter,j) = unitPowerCost(iter,j)*power(iter,j)*tau/sec2hr;  %MWHr
            
            %calculate the amount of data for migration
            
            if((iter<=1) | (s(iter,:) == s(iter-1,:)))        % not equal
                y(j,subscript) = 0;
            else y(j,subscript) = migDataFromDC(iter-1,j);
            end
            m(j,subscript) = m(j,subscript) + w*y(j,subscript);    %the amount of migration data
            if(j == subscript)
                m(j,subscript) = 0;
            end
            r(iter,j) = mu2(iter,j) + beta*sum(m(:,j));
            mm(j,subscript) = unitMigrationCost(j,subscript)*m(j,subscript);
            
            %----Update Queue------
            Q(iter+1,j) = max(Q(iter,j) - x(iter,j),0) + r(iter,j);
            migDataFromDC(iter,j) = min(Q(iter,j), x(iter,j));
        end
        
        sum_migrationCost(a_index,iter) = sum(mm(:,subscript)) + sum(migCostFromFog(iter,:));
        sum_energy(a_index,iter) = sum(energy(iter,:));
        
        %----operation costs------
        cost_matrix(a_index,iter)= sum_bandCost(a_index,iter) + sum_migrationCost(a_index,iter) + sum_energy(a_index,iter)+ sum_energy_fog(a_index,iter);
        
    end  % iter
    
    result_cost(a_index,1) = sum(cost_matrix(a_index,:))/Iters;
    avg_speed(a_index,1) = sum(sum(speed,2))/(Iters*K);
    avg_num(a_index,1) = sum(sum(num,2))/(Iters*K);
    %% ---------------------Length of Queues-------------------------------------
    for iter2 = 1:Iters
        G_avg(a_index,iter2) = sum(G(iter2,:))/M;
        Z_avg(a_index,iter2) = sum(Z(iter2,:))/M;
        Q_avg(a_index,iter2) = sum(Q(iter2,:))/K;
        
        QueLen(a_index,iter2) = (G_avg(a_index,iter2) + Z_avg(a_index,iter2) + Q_avg(a_index,iter2));
    end
    cost_thr(a_index,1) = sum(cost_matrix(a_index,:));
    throughput(a_index,1) = (sum(sum(migDataFromDC)) + sum(sum(fog_B)))/Iters;
    QueLen_bar(a_index,1) = sum(QueLen(a_index,:))/Iters;
end % varialbe V

%% ---------------------Figures-------------------------------------
figure;
set(0,'DefaultTextFontName','Times',...
    'DefaultTextFontSize',16,...
    'DefaultAxesFontName','Times',...
    'DefaultAxesFontSize',12,...
    'DefaultLineLineWidth',1,...
    'DefaultLineMarkerSize',6);
[hAx,hLine1,hLine2] =plotyy(a_options,QueLen_bar(:,1),a_options,result_cost(:,1));
%[hAx,hLine1,hLine2] =plotyy((min_V+step:step:max_V),QueLen_bar(2:end,1),(min_V+step:step:max_V),result_cost(2:end,1));
set(hLine1,'LineStyle','-','Marker','o','Color','b');
set(hLine2,'LineStyle','-','Marker','s','Color','r');
set(hAx(1),'YColor','b');
set(hAx(2),'YColor','r');
grid
%set(hAx(1),'ylim',[15,23]);
% set(hAx(2),'ylim',[10,11.4]);
% %set(hAx(1),'yTick',15:2:23);
% set(hAx(2),'yTick',10:0.2:11.4);
xlabel('V');
ylabel(hAx(1),'Time Aavraged Queue Size (GB)');
ylabel(hAx(2),'Time Aavraged Operation Cost ($)');
hl = legend('Queue Size','Cost');
set(hl,'Location','Best');
