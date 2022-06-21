% Function to dealise a Fourier transformed vector
% Zheng Zheng, June 2022

% input: a Fourier transformed vector
% output: dealised vector
function [dealised_vector] = dealising(vector)
    l = length(vector);
    mode = floor(l/2); % middle point 
    kill = floor(l/6);
    vector(mode - kill : mode + kill) = 0;
    dealised_vector = vector;
end


