% 
% Zheng Zheng, June 2022

% input: 
% output:  
function [x_new, y_new] = inverse_matrix(a, b, c, d, e, f)
    A = [a, b ; c, d];
    B = [e;f];
    X = A\B;
    x_new = X(1);
    y_new = X(2);
end


