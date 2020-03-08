function [Breg_Err, Fvalue_Err] = evaluate(X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code output two optimality measures for an iterate.
% ---------input----------
% X: any stacked iterate;
% ---------output---------
% Breg_Err: the Bregman distance optimality error;
% Fvalue_Err: the function value optimality error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global U V Num_Nodes F_opt X_opt
%%% Bregman distance optimality gap %%%%
Breg_Err = 0;
for i = 1:Num_Nodes
    Breg_Err = Breg_Err + F(U{i}, V{i}, X(i,:)) ...
        - grad(U{i}, V{i}, X_opt) * X(i,:)';
end
Breg_Err = abs(Breg_Err - F_opt);


%%% function value ptimality gap %%%%
U_stack = cat(1, U{:});
V_stack = cat(1, V{:});
Fvalue_Err = 0;
for i = 1:Num_Nodes
    Fvalue_Err = Fvalue_Err + F(U_stack, V_stack, X(i,:)) - F_opt;
end
Fvalue_Err = abs(Fvalue_Err)/Num_Nodes;