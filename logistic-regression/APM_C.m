function [Err, cost_counter] = APM_C(Lap, L, mu, beta0, is_strongly_convex, inner_consensus_const)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code describes the algorithm--APM_C (Algorithm 2 in [1]).
% ---------input----------
% Lap:     a Laplacian matrix;
% L:       Lipschitz constant;
% mu:      strong convexity constant which is not required if 
%           is_strongly_convex == false;
% beta0:   constant;
% is_strongly_convex: logical value;
% inner_consensus_const: constant for determining T_k;
% ---------output---------
% Err: a two-row matrix with the first row recoding the Bregman distance 
%      optimality gap and the second row the function value optimality gap.
% cost_counter: a three-row matrix with frst row recording the total cost; 
%               the second row communication cost and the third row
%               gradient computation cost.
% --------reference-------
% [1] Li, Huan, Cong Fang, Wotao Yin, and Zhouchen Lin. "A Sharp Convergence
%     Rate Analysis for Distributed Accelerated Gradient Methods." arXiv 
%     preprint arXiv:1810.01053 (2018).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global U V col X0 Num_Nodes Niter comp_time_unit comm_time_unit

Err = zeros(Niter+1, 2);

%%%% Weight matrix %%%%
Lap_eig    = sort(eig(Lap));
W          = eye(Num_Nodes) - Lap/Lap_eig(end);
W_eig      = sort(eig(W));

%%%% Parameter %%%%
eta = (1 - sqrt(1-W_eig(end-1)^2))/(1 + sqrt(1-W_eig(end-1)^2));
vartheta               = zeros(1, Niter);
inner_consensus_factor = zeros(1, Niter);
if ~is_strongly_convex
    mu = 0;  theta = zeros(1, Niter+1);  theta(1) = 1;
    for k = 1:Niter
        theta(k+1)  = (sqrt(theta(k)^4 + 4*theta(k)^2) - theta(k)^2)/2;
        vartheta(k) = theta(k+1)^2;
        inner_consensus_factor(k) = log(k) / sqrt(1-W_eig(end-1));
    end
else
    theta = ones(1, Niter+1) * sqrt(mu/L);
    for k = 1:Niter
        vartheta(k)               = (1 - sqrt(mu/L))^(k+1);
        inner_consensus_factor(k) = k * sqrt(mu/L) / sqrt(1-W_eig(end-1));
    end
end

% the first row is for total cost; the second row is for communication
% cost; the third row is for gradient computation cost.
cost_counter      = zeros(3, Niter+1); 

%%%% Update %%%%
old_X = X0;  X = X0;
[Err(1,1), Err(1,2)] = evaluate(X);
for k = 1:Niter
    Y = X + (L*theta(k+1) - mu)/(L - mu) * (1-theta(k))/theta(k) * (X - old_X);
    Grad_Y  = zeros(Num_Nodes, col);
    for i = 1:Num_Nodes
        Grad_Y(i,:) = grad(U{i}, V{i}, Y(i,:));
    end
    initial_Z = Y - Grad_Y/L;
    Z = initial_Z;  old_Z = initial_Z;
    T = ceil(inner_consensus_const * inner_consensus_factor(k));
    for inner_consensus_counter = 1:T
        temp = (1+eta) * W * Z - eta * old_Z;
        old_Z = Z;  Z = temp;
    end
    temp = (L*vartheta(k) * initial_Z + beta0 * Z)/(L*vartheta(k) + beta0);
    old_X = X;  X = temp;
    %%% Evaluate %%%
    [Err(k+1, 1), Err(k+1, 2)] = evaluate(X); 
    cost_counter(1,k+1) = cost_counter(1,k) + T*comm_time_unit + 1*comp_time_unit;
    if mod(k, 100) == 0
        fprintf('APM_C: %d-th iteration, the error is %f\n', k, Err(k+1,1));
    end
end
cost_counter(3,:) = (0:Niter)*1*comp_time_unit;
cost_counter(2,:) = cost_counter(1,:) - cost_counter(3,:);
