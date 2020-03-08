function [Err, cost_counter] = DIGing(Lap, alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code describes the algorithm--DIGing/NEXT (Algorithm 1 in [1] but 
% with an adapt-then-combine version and Algorithm 1 in [2]).
% ---------input----------
% Lap:     A Laplacian matrix;
% alpha:   Step size (which corresponds to I*alpha/tau in [1] if the 
%            surrogate function is chosen as (13))
% ---------output---------
% Err: a two-row matrix with the first row recoding the Bregman distance 
%      optimality gap and the second row the function value optimality gap.
% cost_counter: a three-row matrix with frst row recording the total cost; 
%               the second row communication cost and the third row
%               gradient computation cost.
% --------reference-------
% [1] Nedic, Angelia, Alex Olshevsky, and Wei Shi. "Achieving geometric 
%     convergence for distributed optimization over time-varying graphs." 
%     SIAM Journal on Optimization 27, no. 4 (2017): 2597-2633.
% [2] Di Lorenzo, Paolo, and Gesualdo Scutari. "Next: In-network nonconvex
%     optimization." IEEE Transactions on Signal and Information Processing 
%     over Networks 2, no. 2 (2016): 120-136.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global A b col X0 Niter Num_Nodes comp_time_unit comm_time_unit

Err = zeros(Niter+1, 2);
W = eye(Num_Nodes) - Lap;

%%%% Initialization %%%%
X     = X0;
Y     = zeros(Num_Nodes, col);
for i = 1:Num_Nodes
    Y(i,:) = 2*(X(i,:)*A{i}'-b{i}')*A{i};
end
[Err(1,1), Err(1,2)] = evaluate(X);

for k = 1:Niter
    %%% Update %%%
    X_new = W * (X - alpha * Y);
    Y     = W * Y;
    for i = 1:Num_Nodes
        Y(i,:) = Y(i,:) + 2 * (X_new(i,:) - X(i,:)) * A{i}' * A{i};
    end 
    X     = X_new;   
    [Err(k+1, 1), Err(k+1, 2)] = evaluate(X); 
    if mod(k, 100) == 0
        fprintf('NEXT: %d-th iteration, the error is %f\n', k, Err(k+1,1));
    end
end
% the first row is for total cost; the second row is for communication
% cost; the third row is for gradient computation cost.
cost_counter      = zeros(3, Niter+1);
cost_counter(2,:) = (0:Niter)*2*comm_time_unit;
cost_counter(3,:) = (0:Niter)*1*comp_time_unit;
cost_counter(1,:) = cost_counter(2,:) + cost_counter(3,:);