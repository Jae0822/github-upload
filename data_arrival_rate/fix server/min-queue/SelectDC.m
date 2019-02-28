function [index, output]=SelectDC(DC_que,migration_data,migration_cost,V)
minVal = +inf;
beta  = 0.01;
num_DC = size(DC_que,2);
output = zeros(1,num_DC);
for j=1:num_DC
    temp = beta*DC_que(j)*sum(migration_data(:,j)) + V*sum(migration_cost(:,j));
    if(temp<minVal)
        minVal = temp;
        index = j;
    end
end
output(1,index) = 1;
end