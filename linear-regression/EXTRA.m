function [Err, cost_counter] = EXTRA(Lap, alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code describes the algorithm--EXTRA (Algorithm 1 in [1]).
% ---------input----------
% Lap:     a Laplacian matrix;
% alpha:   step size
% ---------output---------
% Err: a two-row matrix with the first row recording the Bregman distance 
%      optimality gap and the second row the function value optimality gap.
% cost_counter: a three-row matrix with frst row recording the total cost; 
%               the second row communication cost and the third row
%               gradient computation cost.
% --------reference-------
% [1] Shi, Wei, Qing Ling, Gang Wu, and Wotao Yin. "Extra: An exact 
%     first-order algorithm for decentralized consensus optimization." SIAM 
%     Journal on Optimization 25, no. 2 (2015): 944-966.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global A b col X0 Niter Num_Nodes comp_time_unit comm_time_unit

Err = zeros(Niter+1, 2);
W = eye(Num_Nodes) - Lap;
W_tilde = (eye(Num_Nodes) + W)/2;

%%%% Initialization %%%%
X     = X0;
[Err(1,1), Err(1,2)] = evaluate(X);
Grad  = zeros(Num_Nodes, col);
for i = 1:Num_Nodes
    Grad(i,:) = 2*(X(i,:)*A{i}'-b{i}')*A{i};
end
%%%% first step update %%%%
X_new    = W * X - alpha * Grad; 
[Err(2,1), Err(2,2)] = evaluate(X_new);

for k = 2:Niter
    %%% Update %%%
    Grad_new = zeros(Num_Nodes, col);
    for i = 1:Num_Nodes
        Grad_new(i,:) = 2*(X_new(i,:)*A{i}'-b{i}')*A{i};
    end
    temp = (eye(Num_Nodes) + W) * X_new - W_tilde * X ...
        - alpha * (Grad_new - Grad);
    X = X_new;  X_new = temp;  Grad = Grad_new;
    [Err(k+1, 1), Err(k+1, 2)] = evaluate(X_new); 
    if mod(k, 100) == 0
        fprintf('EXTRA: %d-th iteration, the error is %f\n', k, Err(k+1,1));
    end
end
% the first row is for total cost; the second row is for communication
% cost; the third row is for gradient computation cost.
cost_counter          = zeros(3, Niter+1);
cost_counter(2,2:end) = ((1:Niter)*2 - 1) * comm_time_unit;
cost_counter(3,2:end) = (1:Niter)*1*comp_time_unit;
cost_counter(1,:)     = cost_counter(2,:) + cost_counter(3,:);