% Recurrent "flow" analysis 
% Zheng Zheng, June 2022

function [diff] = recurrent_flow(tmax,dt,xyz)
    diff = zeros(floor(tmax/dt)); 
    for time = 1 : (floor(tmax/dt))
        for i = 1 : time    
            res = sqrt((xyz(time,1) - xyz(i,1))^2  + (xyz(time,2) - xyz(i,2))^2 + (xyz(time,3) - xyz(i,3))^2);
            diff(i, time) = res;
        end
    end
end


