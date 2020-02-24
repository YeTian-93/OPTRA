function g = grad(U, V, x) 
global r

total = size(U,1);
R = total * r;
g = 0;
for m = 1:total
    g = g - V(m)/(1 + exp(V(m)*U(m,:)*x')) * U(m,:);
end
g = g + 2*R*x; 
end