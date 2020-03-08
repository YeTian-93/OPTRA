function g = stoc_grad(i, x, batch_size_portion) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code outputs a stochastic gradient corresponding to a random
% mini-batch sampled uniformly from the local data set.
% ---------input----------
% i:                  index of the agent;
% x:                  current iterate;
% batch_size_portion: the propotion of the size of the mini-batch to the
%                     full local data set.
% ---------output---------
% g:                  resulted stochastic gradient.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global A b

feature = A{i};      label = b{i};
total = size(feature,1);
batch_size = ceil(total * batch_size_portion);

rand_sp = randsample(total, batch_size);
g = 0;
for k = 1:batch_size
    m = rand_sp(k);    
    g = g + 2 * (x * feature(m, :)' - label(m)) * feature(m, :);    
end
g = g / batch_size * total;
end