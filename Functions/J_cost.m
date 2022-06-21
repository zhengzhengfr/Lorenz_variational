% Function to calculate the scalar residual J
% Zheng Zheng, June 2022

% input: residual field vectors in x, y and z
% output: scalar J
function [J] = J_cost(r1, r2, r3)
    % Integration from 0 to 1 
    res1 = mean(r1.^2); 
    res2 = mean(r2.^2); 
    res3 = mean(r3.^2); 
    J = sqrt(res1 + res2 + res3);
end


