function [Err, cost_counter] = DPSGD(Lap, step, batch_size_portion)
global col X0 Niter Num_Nodes 

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
cost_counter = zeros(3, num_iter+1);
cost_counter(1,:) = (0:num_iter)*(1+batch_size_portion);
cost_counter(2,:) = (0:num_iter)*1;
cost_counter(3,:) = (0:num_iter)*batch_size_portion;