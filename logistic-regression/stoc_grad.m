function g = stoc_grad(i, x, batch_size_portion) 
global U V

feature = U{i};      label = V{i};
total = size(feature,1);
batch_size = ceil(total * batch_size_portion);

rand_sp = randsample(total, batch_size);
g = 0;
for k = 1:batch_size
    m = rand_sp(k);    
    g = g - label(m)/(1 + exp(label(m)*feature(m,:)*x')) * feature(m,:);
end
g = g / batch_size * total;
end