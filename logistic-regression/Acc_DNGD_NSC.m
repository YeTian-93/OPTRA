function [Err, cost_counter] = Acc_DNGD_NSC(Lap, L, eta, step_size_mode, t0, beta)
global U V col Niter Num_Nodes

Err = zeros(Niter+1, 2);
W   = eye(Num_Nodes) - Lap;

%%%% Initialization %%%%
X     = zeros(Num_Nodes, col);
Z     = zeros(Num_Nodes, col);  % Z is to denote V in the pseudocode.
Y     = zeros(Num_Nodes, col);
S     = zeros(Num_Nodes, col);

for i = 1:Num_Nodes
    S(i,:) = grad(U{i}, V{i}, X(i,:));
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
    Z = W * Z - Eta(k)/alpha(k) * S;
    new_Y = (1-alpha(k+1)) * X + alpha(k+1) * Z;
    S     = W * S;
    for i = 1:Num_Nodes
        S(i,:) = S(i,:) + grad(U{i}, V{i}, new_Y(i,:))...
            - grad(U{i}, V{i}, Y(i,:));
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
cost_counter(1,:) = (0:Niter)*4;
cost_counter(2,:) = (0:Niter)*3;
cost_counter(3,:) = (0:Niter)*1;