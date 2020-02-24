function f = F(U, V, x) 
global r 

total = size(U,1);
R = r * total;
f = 0;
for m = 1:total
    f = f + log(1 + exp(-V(m)*U(m,:)*x'));
end
f = f + R * norm(x)^2;
end