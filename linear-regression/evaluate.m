function [Breg_Err, Fvalue_Err] = evaluate(X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code output two optimality measures for an iterate.
% ---------input----------
% X: any stacked iterate;
% ---------output---------
% Breg_Err: the Bregman distance optimality error;
% Fvalue_Err: the function value optimality error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global A b X_opt fmin Num_Nodes

%%%% Bregman distance optimality gap %%%%
Breg_Err = 0;
for i = 1:Num_Nodes
    Breg_Err = Breg_Err + norm(A{i} * X(i,:)' - b{i})^2 ...
        - X(i,:) * 2 * A{i}' * (A{i}*X_opt-b{i});
end
Breg_Err = abs(Breg_Err - fmin);


%%%% function value ptimality gap %%%%
A_stack = cat(1, A{:});
b_stack = cat(1, b{:});
Fvalue_Err = 0;
for i = 1:Num_Nodes
    Fvalue_Err = Fvalue_Err + norm(A_stack * X(i,:)' - b_stack)^2 - fmin;
end
Fvalue_Err = abs(Fvalue_Err)/Num_Nodes; 