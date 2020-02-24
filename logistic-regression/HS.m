function h = HS(U, V, x) 

[total, col] = size(U);
H = zeros(col);
for m = 1:total
    temp = exp(V(m)*U(m,:)*x');
    H = H + temp/(1+temp)^2 * V(m)^2 * U(m,:)' * U(m,:);
end
h = eigs(H, 1);
end