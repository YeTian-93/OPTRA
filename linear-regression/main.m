clear all;
close all;
rng(052019,'v4');

global A b fmin X_opt col row Num_Nodes Niter X0

%%%% Problem parameters %%%%
Num_Nodes  = 20;
row        = 10;
col        = 500;
Niter      = 20000;

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

%% OPTRA_C
nu = 100;
Num_InnerConsensus = 4;
% This parameter denotes the number of communications between two
% consecutive gradient evalutions in OptPrimalDual and should be an
% integer not less than 2.
[OPTRA_err, OPTRA_counter, new_eigengap] = OPTRA_C(Lap, L_f, nu, Num_InnerConsensus);
disp(['our method: ', num2str(OPTRA_err(end,1))]);

% %% theoretical upper bound
% R_x = norm(X0 - ones(Num_Nodes,1)*X_opt', 'fro')^2;  
% R_y = 0;
% for i = 1:Num_Nodes
%     R_y = R_y + norm(2*A{i}'*(A{i}*X_opt-b{i}))^2;
% end
% upperbound = 2*L_f*R_x./(1:Niter+1).^2 + ...
%     (2/nu*R_x + 2*nu*R_y/new_eigengap)./(1:Niter+1);
fprintf('The new_eigengap is %f.\n', new_eigengap);


%%
% ---- the definition of the input arguments and outputs of all the above
% ---- algorithms can be found in the corresponding .m files
toc

% DIGing_err              = zeros(Niter+1, 2);
% EXTRA_err               = zeros(Niter+1, 2);
% Acc_DNGD_NSC_err        = zeros(Niter+1, 2);
% DIGing_counter          = zeros(3, Niter+1);
% EXTRA_counter           = zeros(3, Niter+1);
% Acc_DNGD_NSC_counter    = zeros(3, Niter+1);

% save result.mat

load result.mat

%% plot figures
figure
subplot(1,3,1)
fmean = semilogy(DIGing_counter(1,:), DIGing_err(:,1), EXTRA_counter(1,:), EXTRA_err(:,1),...
    Acc_DNGD_NSC_counter(1,:), Acc_DNGD_NSC_err(:,1),...
    APM_C_counter(1,:), APM_C_err(:,1), DPSGD_counter(1,:), DPSGD_err(:,1), ...
    OPTRA_counter(1,:), OPTRA_err(:,1));
hold on
set(fmean, 'linewidth', 2);
xlabel({'total cost'}, 'FontSize',16)
ylabel({'Bregman distance optimality gap'}, 'FontSize',16)
legend({'DIGing/NEXT', 'EXTRA', 'Acc-DNGD-NSC', 'APM-C', 'DPSGD', 'Our method'}, 'FontSize',16)
xlim([0 DPSGD_counter(1,end)])
set(gca,'FontSize',16)

subplot(1,3,2)
fmean = semilogy(DIGing_counter(2,:), DIGing_err(:,1), EXTRA_counter(2,:), EXTRA_err(:,1),...
    Acc_DNGD_NSC_counter(2,:), Acc_DNGD_NSC_err(:,1),...
    APM_C_counter(2,:), APM_C_err(:,1), DPSGD_counter(2,:), DPSGD_err(:,1), ...
    OPTRA_counter(2,:), OPTRA_err(:,1));
hold on
set(fmean, 'linewidth', 2);
xlabel({'communication cost'}, 'FontSize',16)
% ylabel({'Bregman distance optimality gap'}, 'FontSize',20)
legend({'DIGing/NEXT', 'EXTRA', 'Acc-DNGD-NSC', 'APM-C', 'DPSGD', 'Our method'}, 'FontSize',16)
xlim([0 DIGing_counter(2,end)])
set(gca,'FontSize',16)

subplot(1,3,3)
fmean = semilogy(DIGing_counter(3,:), DIGing_err(:,1), EXTRA_counter(3,:), EXTRA_err(:,1),...
    Acc_DNGD_NSC_counter(3,:), Acc_DNGD_NSC_err(:,1),...
    APM_C_counter(3,:), APM_C_err(:,1), DPSGD_counter(3,:), DPSGD_err(:,1), ...
    OPTRA_counter(3,:), OPTRA_err(:,1));
hold on
% semilogy(OPTRA_counter(3,:), upperbound, 'ro');
set(fmean, 'linewidth', 2);
xlabel({'computation cost'}, 'FontSize',16)
% ylabel({'Bregman distance optimality gap'}, 'FontSize',20)
legend({'DIGing/NEXT', 'EXTRA', 'Acc-DNGD-NSC', 'APM-C', 'DPSGD', 'Our method'}, 'FontSize',16)
xlim([0 DPSGD_counter(3,end)])
ylim([1e-8 inf])
set(gca,'FontSize',16)


figure
subplot(1,3,1)
fmean = semilogy(DIGing_counter(1,:), DIGing_err(:,2), EXTRA_counter(1,:), EXTRA_err(:,2),...
    Acc_DNGD_NSC_counter(1,:), Acc_DNGD_NSC_err(:,2),...
    APM_C_counter(1,:), APM_C_err(:,2), DPSGD_counter(1,:), DPSGD_err(:,2), ...
    OPTRA_counter(1,:), OPTRA_err(:,2));
hold on
set(fmean, 'linewidth', 2);
xlabel({'total cost'}, 'FontSize',16)
ylabel({'function value optimality gap'}, 'FontSize',16)
legend({'DIGing/NEXT', 'EXTRA', 'Acc-DNGD-NSC', 'APM-C', 'DPSGD', 'Our method'}, 'FontSize',16)
xlim([0 DPSGD_counter(1,end)])
set(gca,'FontSize',16)

subplot(1,3,2)
fmean = semilogy(DIGing_counter(2,:), DIGing_err(:,2), EXTRA_counter(2,:), EXTRA_err(:,2),...
    Acc_DNGD_NSC_counter(2,:), Acc_DNGD_NSC_err(:,2),...
    APM_C_counter(2,:), APM_C_err(:,2), DPSGD_counter(2,:), DPSGD_err(:,2), ...
    OPTRA_counter(2,:), OPTRA_err(:,2));
hold on
set(fmean, 'linewidth', 2);
xlabel({'communication cost'}, 'FontSize',16)
% ylabel({'function value optimality gap'}, 'FontSize',20)
legend({'DIGing/NEXT', 'EXTRA', 'Acc-DNGD-NSC', 'APM-C', 'DPSGD', 'Our method'}, 'FontSize',16)
xlim([0 DIGing_counter(2,end)])
set(gca,'FontSize',16)

subplot(1,3,3)
fmean = semilogy(DIGing_counter(3,:), DIGing_err(:,2), EXTRA_counter(3,:), EXTRA_err(:,2),...
    Acc_DNGD_NSC_counter(3,:), Acc_DNGD_NSC_err(:,2),...
    APM_C_counter(3,:), APM_C_err(:,2), DPSGD_counter(3,:), DPSGD_err(:,2), ...
    OPTRA_counter(3,:), OPTRA_err(:,2));
hold on
set(fmean, 'linewidth', 2);
xlabel({'computation cost'}, 'FontSize',16)
% ylabel({'function value optimality gap'}, 'FontSize',20)
legend({'DIGing/NEXT', 'EXTRA', 'Acc-DNGD-NSC', 'APM-C', 'DPSGD', 'Our method'}, 'FontSize',16)
xlim([0 DPSGD_counter(3,end)])
set(gca,'FontSize',16)