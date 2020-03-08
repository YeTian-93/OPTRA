function f = F(U, V, x) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code outputs the function value corresponding to the data set (U, V).
% ---------input----------
% U:                  feature matrix;
% V:                  label vector;
% x:                  current iterate.
% ---------output---------
% f:                  resulted function value.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global r 

total = size(U,1);
R = r * total;
f = 0;
for m = 1:total
    f = f + log(1 + exp(-V(m)*U(m,:)*x'));
end
f = f + R * norm(x)^2;
end