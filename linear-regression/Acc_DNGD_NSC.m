function [Err, cost_counter] = Acc_DNGD_NSC(Lap, L, eta, step_size_mode, t0, beta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code describes the algorithm--Acc-DNGD-NSC in [1].
% ---------input----------
% Lap:     a Laplacian matrix;
% L:       Lipschitz constant;
% eta:     step size;
% step_size_mode: the default mode is constant step size.  The user need to 
%                 input "dimi" for this argument to invoke diminishing 
%                 step size mode;
% t0, beta: constants needed in dimininishing step size and is not required
%              for constant step size mode
% ---------output---------
% Err: a two-row matrix with the first row recoding the Bregman distance 
%      optimality gap and the second row the function value optimality gap.
% cost_counter: a three-row matrix with frst row recording the total cost; 
%               the second row communication cost and the third row
%               gradient computation cost.
% --------reference-------
% [1] Qu, Guannan, and Na Li. "Accelerated distributed Nesterov gradient 
%     descent." arXiv preprint arXiv:1705.07176 (2017).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global A b col Niter Num_Nodes comp_time_unit comm_time_unit

Err = zeros(Niter+1, 2);
W   = eye(Num_Nodes) - Lap;

%%%% Initialization %%%%
X     = zeros(Num_Nodes, col);
V     = zeros(Num_Nodes, col);
Y     = zeros(Num_Nodes, col);
S     = zeros(Num_Nodes, col);

for i = 1:Num_Nodes
    S(i,:) = 2*(X(i,:)*A{i}'-b{i}')*A{i};
end
[Err(1,1), Err(1,2)] = evaluate(Y);

Eta = eta * ones(1, Niter+1);
switch step_size_mode
    case 'dimi'
        for k = 1:Niter+1
            Eta(k) = eta/(k - 1 + t0)^beta;
        end
end

alpha = zeros(1, Niter+1);
alpha(1) = sqrt(Eta(1)*L);

for k = 1:Niter
    temp = Eta(k+1)/Eta(k)*alpha(k)^2;
    alpha(k+1) = (sqrt(temp^2 + 4*temp)-temp)/2;
    %%% Update %%%
    X = W * Y - Eta(k) * S;
    V = W * V - Eta(k)/alpha(k) * S;
    new_Y = (1-alpha(k+1)) * X + alpha(k+1) * V;
    S     = W * S;
    for i = 1:Num_Nodes
        S(i,:) = S(i,:) + 2 * (new_Y(i,:) - Y(i,:)) * A{i}' * A{i};
    end 
    Y     = new_Y;
    
    [Err(k+1, 1), Err(k+1, 2)] = evaluate(Y);  
    if mod(k, 100) == 0
        fprintf('Acc_DNGD: %d-th iteration, the error is %f\n', k, Err(k+1,1));
    end
end
% the first row is for total cost; the second row is for communication
% cost; the third row is for gradient computation cost.
cost_counter      = zeros(3, Niter+1); 
cost_counter(2,:) = (0:Niter)*3*comm_time_unit;
cost_counter(3,:) = (0:Niter)*1*comp_time_unit;
cost_counter(1,:) = cost_counter(2,:) + cost_counter(3,:);