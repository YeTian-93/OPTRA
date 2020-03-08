function data_generator(N, step_size)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code generates data for the decentralized logistic regression 
% problem on the Parkinson?s Disease Classification Data Set (available at
% https://archive.ics.uci.edu/ml/datasets/Parkinson%27s+Disease+Classification).
% See section F.2 in our paper for details of the generation process.
% ---------input----------
% N:          number of agents;
% step_size:  the step size used in the centralized gradient descent
%             algorithm for finding an x^star and the optimal function value.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(112019,'v4');
global Num_Nodes 

%%%%% preprocess data %%%%%
load('data_pd.mat')
data                  = pdspeechfeatures{:,:};
features              = data(:, 2:end-1);
[num_sp, col]         = size(features);
V_stack               = data(:,end);
U_stack               = (features-ones(num_sp,1)*min(features))...
                         ./(ones(num_sp,1)*(max(features)-min(features)));
V_stack(V_stack == 0) = -1;

%%%%% distribute data evenly %%%%%
U = cell(1, Num_Nodes);
V = cell(1, Num_Nodes);
row = floor(num_sp/Num_Nodes);
for i =1:Num_Nodes-1
    U{i} = U_stack((i-1)*row+1:i*row, :);
    V{i} = V_stack((i-1)*row+1:i*row);
end
U{Num_Nodes} = U_stack((Num_Nodes-1)*row+1:end, :);
V{Num_Nodes} = V_stack((Num_Nodes-1)*row+1:end, :);

%%%%% minimizing by centralized gradient descent %%%%%
x        = zeros(1, col);
gradnorm = zeros(N,1);
F_value  = zeros(N,1);
F_desc   = zeros(N-1,1);

for k = 1:N
    g           = grad(U_stack, V_stack, x);
    x           = x - step_size * g;
    gradnorm(k) = norm(g);
    F_value(k)  = F(U_stack, V_stack, x);
    if k > 1
        F_desc(k-1) = abs(F_value(k-1) - F_value(k));
    end
    if mod(k, 100) == 0
        fprintf('gradient descent: %d-th iteration, the error is %f\n', k, gradnorm(k));
    end
end

F_opt = F_value(N);
X_opt = x;

%%%%% plot %%%%%
figure;

subplot(1,3,1);
semilogy(1:N, gradnorm);
xlabel('Iteration')
ylabel('Gradient Norm')

subplot(1,3,2);
semilogy(1:N, F_value);
xlabel('Iteration')
ylabel('F value')

subplot(1,3,3);
semilogy(1:N-1, F_desc);
xlabel('Iteration')
ylabel('descent per step of F value')

save('data.mat', 'U', 'V', 'F_opt', 'X_opt', 'F_value', 'col', 'row')
fprintf('Data is generated\n');
end