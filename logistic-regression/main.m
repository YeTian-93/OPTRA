clear all;
close all;
rng(052019,'v4');

global U V F_opt X_opt row col Niter X0 Num_Nodes r comp_time_unit comm_time_unit

%%%% Problem parameters %%%%
Num_Nodes      = 60;
r              = 0; % regularization constant
comp_time_unit = 1;
comm_time_unit = 1;

% data_generator(150000, 0.05);
% Network_Laplician(Num_Nodes);
% ---- the user can use the above two lines of code to generate the problem
% ---- data and the graph Laplacian.
load data.mat
load laplacian.mat
Niter         = 15000;

tic
L_f = 24;
X0 = rand(Num_Nodes, col);
%% Optimal Primal Dual
nu = 1500;   
Num_InnerConsensus = 4;  
% This parameter denotes the number of communications between two 
% consecutive gradient evalutions in OptPrimalDual and should be an
% integer not less than 2.
[OptPrimalDual_err, OptPrimalDual_counter] = OPTRA(Lap, L_f, nu, Num_InnerConsensus);
disp(['our method: ', num2str(OptPrimalDual_err(end,1))]);

%% Acc_DNGD_NSC
eta = 0.01/L_f;
[Acc_DNGD_NSC_err, Acc_DNGD_NSC_counter] = Acc_DNGD_NSC(Lap, L_f, eta, 'const', [], []);
disp(['Acc_DNGD_NSC: ', num2str(Acc_DNGD_NSC_err(end,1))]);

%% APM_C
[APM_C_err, APM_C_counter] = APM_C(Lap, L_f, [], 1e4, false, 0.2);
disp(['APM_C: ', num2str(APM_C_err(end,1))]);

%% NEXT
alpha = 0.01;
[NEXT_err, NEXT_counter] = NEXT(Lap, alpha);
disp(['NEXT: ', num2str(NEXT_err(end,1))]);

%% EXTRA
alpha = 0.005;
[EXTRA_err, EXTRA_counter] = EXTRA(Lap, alpha);
disp(['EXTRA: ', num2str(EXTRA_err(end,1))]);

%% DPSGD
step = 0.001;
batch_size_portion = 0.2;
[DPSGD_err, DPSGD_counter] = DPSGD(Lap, step, batch_size_portion);
disp(['DPSGD: ', num2str(DPSGD_err(end,1))]);

%%
% ---- the definition of the input arguments and outputs of all the above
% ---- algorithms can be found in the corresponding .m files
toc

save final-result.mat