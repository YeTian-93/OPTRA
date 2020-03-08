function [Err, cost_counter] = DPSGD(Lap, step, batch_size_portion)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code describes the algorithm--DPSGD in [1].
% ---------input----------
% Lap:                a Laplacian matrix;
% step:               step size;
% batch_size_portion: the propotion of the size of the mini-batch to the
%                     full local data set.
% ---------output---------
% Err: a two-row matrix with the first row recoding the Bregman distance 
%      optimality gap and the second row the function value optimality gap.
% cost_counter: a three-row matrix with frst row recording the total cost; 
%               the second row communication cost and the third row
%               gradient computation cost.
% --------reference-------
% [1] Lian, X., Zhang, C., Zhang, H., Hsieh, C.-J., Zhang, W., and Liu, J. 
% (2017). Can decentralized algorithms outperform centralized algorithms? a
% case study for decentralized parallel stochastic gradient descent. In
% Advances in Neural Information Processing Systems, pages 5330?5340.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global col X0 Niter Num_Nodes comp_time_unit comm_time_unit

num_iter = ceil(2 *Niter);
Err = zeros(num_iter+1, 2);
W = eye(Num_Nodes) - Lap;

%%%% Initialization %%%%
X = X0;
Y = zeros(Num_Nodes, col);

[Err(1,1), Err(1,2)] = evaluate(X);

for k = 1:num_iter
    for i=1:Num_Nodes
        Y(i,:) = stoc_grad(i, X(i,:), batch_size_portion);
    end
    
    X = W * X - step * Y;
    
    
    [Err(k+1, 1), Err(k+1, 2)] = evaluate(X);
    if mod(k, 100) == 0
        fprintf('DPSGD: %d-th iteration, the error is %f\n', k, Err(k+1,1));
    end
    
end
% the first row is for total cost; the second row is for communication
% cost; the third row is for gradient computation cost.
cost_counter      = zeros(3, num_iter+1);
cost_counter(2,:) = (0:num_iter)*1*comm_time_unit;
cost_counter(3,:) = (0:num_iter)*batch_size_portion*comp_time_unit;
cost_counter(1,:) = cost_counter(2,:) + cost_counter(3,:);