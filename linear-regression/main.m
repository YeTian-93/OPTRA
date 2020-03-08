clear all;
close all;
rng(052019,'v4');

global A b fmin X_opt col row Num_Nodes Niter X0 comp_time_unit comm_time_unit

%%%% Problem parameters %%%%
Num_Nodes      = 20;
row            = 10;
col            = 500;
Niter          = 20000;
comp_time_unit = 5;
comm_time_unit = 1;

% LeastSquare_Data_Gen(Num_Nodes, 0.95)
% Network_Laplician(Num_Nodes);
% ---- the user can use the above two lines of code to generate the problem
% ---- data and the graph Laplacian.
load data.mat
load laplacian.mat

if Num_Nodes*row <= col
    condition_num = F_Hessian_eigs(end)/F_Hessian_eigs(col-Num_Nodes*row + 1);
else
    condition_num = F_Hessian_eigs(end)/F_Hessian_eigs(1);
end
fprintf('The condition number of the problem is %f.\n', condition_num);
X0 = rand(Num_Nodes, col);

tic
%% DPSGD
step = 1e-5;
batch_size_portion = 0.2;
[DPSGD_err, DPSGD_counter] = DPSGD(Lap, step, batch_size_portion);
disp(['DPSGD: ', num2str(DPSGD_err(end,1))]);

%% DIGing/NEXT
alpha = 1e-5;
[DIGing_err, DIGing_counter] = DIGing(Lap, alpha);
disp(['DIGing: ', num2str(DIGing_err(end,1))]);

%% EXTRA
alpha = 1e-5;
[EXTRA_err, EXTRA_counter] = EXTRA(Lap, alpha);
disp(['EXTRA: ', num2str(EXTRA_err(end,1))]);

%% Acc_DNGD_NSC
eta = 0.005/L_f;
[Acc_DNGD_NSC_err, Acc_DNGD_NSC_counter] = Acc_DNGD_NSC(Lap, L_f, eta, 'const', [], []);
disp(['Acc_DNGD_NSC: ', num2str(Acc_DNGD_NSC_err(end,1))]);

%% APM_C
[APM_C_err, APM_C_counter] = APM_C(Lap, L_f, [], 1e4, false, 0.2);
disp(['APM_C: ', num2str(APM_C_err(end,1))]);

%% OPTRA
nu = 100;
Num_InnerConsensus = 4;
% This parameter denotes the number of communications between two
% consecutive gradient evalutions in OptPrimalDual and should be an
% integer not less than 2.
[OPTRA_err, OPTRA_counter] = OPTRA(Lap, L_f, nu, Num_InnerConsensus);
disp(['our method: ', num2str(OPTRA_err(end,1))]);


%%
% ---- the definition of the input arguments and outputs of all the above
% ---- algorithms can be found in the corresponding .m files
toc

save final-result.mat