function [Err, cost_counter] = OPTRA(Lap, L, nu, Num_InnerConsensus)

global U V row col X0 Niter Num_Nodes 

Err = zeros(Niter+1, 2); % The first column stores optimality error in 
                         % Bregman distance; the second column stores
                         % optimality error in function value difference.

% %%%% Network %%%%
% Lap_eig   = sort(eig(Lap));
% eig_gap   = Lap_eig(2)/Lap_eig(end);
% c1        = (1-sqrt(eig_gap))/(1+sqrt(eig_gap));
% c2        = (1+eig_gap)/(1-eig_gap);
% c3        = 2/(Lap_eig(2)+Lap_eig(end));
% 
% % % Lap       = Lap * 2 / (eigLap(2) + eigLap(end));
% % Lap       = Lap / eigLap(end);
% % J   = eye(Num_Nodes) - ones(Num_Nodes,1) * ones(1,Num_Nodes) / Num_Nodes;
% %%%% Hyperparameter %%%%
% % nu       = 50;
% 
% K = floor(Num_InnerConsensus/2);
% new_Lap = AccGossip(c3*Lap, K, c2);
% new_Lap_largest_eig = 1 + 2 * c1^K/(1 + c1^(2*K));
% gamma   = nu/(nu*L + Niter) ;
% tau     = 4/(nu * Niter * new_Lap_largest_eig);% 1/gamma;
% 
% %%%% Initialization %%%%
% X      = X0;
% H      = X0;  % H is used to denote U in the algorithm
% Y      = zeros(Num_Nodes, col);
% hat_Y  = tau * new_Lap * X;
% 
% [Err(1,1), Err(1,2)] = evaluate(X);

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
theta = 1;
for k = 1:Niter
    new_theta  = (sqrt(theta^4 + 4*theta^2) - theta^2)/2;
    %%% step sizes %%%
    alpha = new_theta/theta - new_theta;
    sigma = 1/new_theta;
    Tau   = tau / theta;
    beta  = theta/new_theta;
    theta = new_theta;   
%     beta  = new_theta/theta;
    
%     theta = new_theta;  
% for k = 1:Niter   
%     %%% step sizes %%%
%     alpha = (k-1)/(k+2);
%     sigma = k/2 + 1;
%     Tau   = tau * (k+1)/2;
%     beta  = (k+2)/(k+1);
    
    %%% Update %%%    
    
    Grad = zeros(Num_Nodes, col);
    for i = 1:Num_Nodes
        Grad(i,:) = grad(U{i}, V{i}, X(i,:));
    end
%     new_H = (eye(Num_Nodes) - new_Lap/new_Lap_largest_eig) ...
%         * (X - gamma * Grad - gamma * hat_Y);
    new_H = (eye(Num_Nodes) - c2 * new_Lap) ...
        * (X - gamma * Grad - gamma * hat_Y);
    X = new_H + alpha * (new_H - H);
    hatX = sigma * X + (1-sigma) * new_H;    
    new_Y = Y + Tau * new_Lap * hatX;
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