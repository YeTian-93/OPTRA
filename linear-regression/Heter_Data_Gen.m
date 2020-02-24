function Heter_Data_Gen(Num_Nodes, omega)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code generates data for the least-square problem.  The data
% generation method is introduced in [1] and summarized in the simulation
% section of our paper.
% ---------input----------
% Num_Nodes:      number of agents;
% omega:          the parameter controlling the condition number of the
%                 problem
% --------reference-------
% [1] Agarwal, Alekh, Sahand Negahban, and Martin J. Wainwright. "Fast 
%     global convergence rates of gradient methods for high-dimensional 
%     statistical recovery." In Advances in Neural Information Processing 
%     Systems, pp. 37-45. 2010.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global row col

Z = normrnd(0, 1, [row*Num_Nodes, col]);
A_stack = zeros(row*Num_Nodes, col);
A_stack(:,1) = Z(:,1)/sqrt(1-omega^2);
for i = 2:col
    A_stack(:,i) = omega * A_stack(:,i-1) + Z(:,i);
end
b_stack       = randn(row*Num_Nodes,1);

A = cell(Num_Nodes, 1);
b = cell(Num_Nodes, 1);
f_Hessian_eigs   = [];

for i = 1:Num_Nodes
    A_stack((i-1)*row+1:i*row,:) = i * A_stack((i-1)*row+1:i*row,:);
    A{i} = A_stack((i-1)*row+1:i*row,:);
    b{i} = b_stack((i-1)*row+1:i*row);
    f_Hessian_eigs  = [f_Hessian_eigs; 2*eig(A{i}'*A{i})];
end

X_opt = A_stack\b_stack;
fmin  = norm(A_stack * X_opt - b_stack)^2;
F_Hessian_eigs   = 2 * sort(eig(A_stack' * A_stack));

L_f = max(f_Hessian_eigs);  % the Lipschitz constant of the augmented function
fprintf('The least-square data is generated.\n')
save('data.mat', 'A', 'b', 'X_opt', 'fmin','F_Hessian_eigs', 'L_f')
end