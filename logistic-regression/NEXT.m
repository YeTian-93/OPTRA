function [Err, cost_counter] = NEXT(Lap, alpha)

global U V col X0 Niter Num_Nodes 

Err = zeros(Niter+1, 2);
W = eye(Num_Nodes) - Lap;

%%%% Initialization %%%%
X = X0;
Y = zeros(Num_Nodes, col);
for i=1:Num_Nodes     
    Y(i,:) = grad(U{i}, V{i}, X(i,:));  
end
[Err(1,1), Err(1,2)] = evaluate(X);

for k = 1:Niter
    
    new_X = W * (X - alpha * Y);
    Y     = W * Y;
    for i = 1:Num_Nodes    
        Y(i,:) = Y(i,:) + grad(U{i}, V{i}, new_X(i,:))...
            - grad(U{i}, V{i}, X(i,:));         
    end
    X = new_X;
    [Err(k+1, 1), Err(k+1, 2)] = evaluate(X); 
    if mod(k, 100) == 0
        fprintf('NEXT: %d-th iteration, the error is %f\n', k, Err(k+1,1));
    end
    
end
% the first row is for total cost; the second row is for communication
% cost; the third row is for gradient computation cost.
cost_counter = zeros(3, Niter+1); 
cost_counter(1,:) = (0:Niter)*3;
cost_counter(2,:) = (0:Niter)*2;
cost_counter(3,:) = (0:Niter)*1;