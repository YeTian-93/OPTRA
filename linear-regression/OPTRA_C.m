function [Err, cost_counter, new_eigengap] = OPTRA_C(Lap, L, nu, Num_InnerConsensus)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code our algorithm--OPTRA-C.
% ---------input----------
% Lap:     a Laplacian matrix;
% L:       lipschitz constant;
% nu:      a constant;
% Num_InnerConsensus:  Number of total inner consensus at each iteration.
% ---------output---------
% Err: a two-row matrix with the first row recoding the Bregman distance 
%      optimality gap and the second row the function value optimality gap.
% cost_counter: a three-row matrix with frst row recording the total cost; 
%               the second row communication cost and the third row
%               gradient computation cost.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global A b col X0 Num_Nodes Niter

Err = zeros(Niter+1, 2); % The first column stores optimality error in 
                         % Bregman distance; the second column stores
                         % optimality error in function value difference.

%%%% Network %%%%
Lap_eig   = sort(eig(Lap));
eig_gap   = Lap_eig(2)/Lap_eig(end);
c0        = (1-sqrt(eig_gap))/(1+sqrt(eig_gap));
c1        = (1+eig_gap)/(1-eig_gap);
c3        = 2/(Lap_eig(2)+Lap_eig(end));

%%%% Hyperparameter %%%%
K         = floor(Num_InnerConsensus/2);
new_Lap   = AccGossip(c3*Lap, K, c1);
c2        = 1/(1 + 2 * c0^K/(1 + c0^(2*K)));
gamma     = nu/(nu*L + Niter+1);
tau       = c2/(nu * (Niter+1));
new_eigengap = (1 - 2 * c0^K/(1 + c0^(2*K)))/(1 + 2 * c0^K/(1 + c0^(2*K)));

%%%% Initialization %%%%
X         = X0;
H         = X0;  % H is to denote U in the pseudocode
Y         = zeros(Num_Nodes, col);
hat_Y     = tau * new_Lap * X;
[Err(1,1), Err(1,2)] = evaluate(X);

%%%% Update %%%%
theta = 1;
for k = 1:Niter
    new_theta  = (sqrt(theta^4 + 4*theta^2) - theta^2)/2;
    %%% step sizes %%%
    alpha = new_theta/theta - new_theta;
    sigma = 1/new_theta;
    Tau   = tau / theta;
    beta  = theta/new_theta;
%     beta  = new_theta/theta;
    
    theta = new_theta;    
    
    Grad  = zeros(Num_Nodes, col);
    for i = 1:Num_Nodes
        Grad(i,:) = 2*(X(i,:)*A{i}'-b{i}')*A{i};        
    end    
    new_H = (eye(Num_Nodes) - c2 * new_Lap) * ...
        (X - gamma * (Grad + hat_Y)); 
%     new_H = X - gamma * (Grad + hat_Y);
    X = new_H + alpha * (new_H - H);
    hat_X = sigma * X + (1-sigma) * new_H;    
    new_Y = Y + Tau * new_Lap * hat_X;
    hat_Y = new_Y + beta * (new_Y - Y);
    H = new_H;  Y = new_Y;   
      
    [Err(k+1, 1), Err(k+1, 2)] = evaluate(X); 
    if mod(k, 100) == 0
        fprintf('OPTRA: %d-th iteration, the error is %f\n', k, Err(k+1,1));
    end
end

% the first row is for total cost; the second row is for communication
% cost; the third row is for gradient computation cost.
cost_counter      = zeros(3, Niter+1); 
cost_counter(1,:) = (0:Niter)*(2*K + 1);
cost_counter(2,:) = (0:Niter)*2*K;
cost_counter(3,:) = (0:Niter)*1;