% NEXT_err                = zeros(Niter+1, 2);
% DIGing_err              = zeros(Niter+1, 2);
% EXTRA_err               = zeros(Niter+1, 2);
% DPSGD_err               = zeros(Niter+1, 2);
% Acc_DNGD_NSC_err        = zeros(Niter+1, 2);
% APM_C_err               = zeros(Niter+1, 2);
% OptPrimalDual_err       = zeros(Niter+1, 2);
% NEXT_counter            = zeros(3, Niter+1);
% DIGing_counter          = zeros(3, Niter+1);
% EXTRA_counter           = zeros(3, Niter+1);
% DPSGD_counter           = zeros(3, Niter+1);
% Acc_DNGD_NSC_counter    = zeros(3, Niter+1);
% APM_C_counter           = zeros(3, Niter+1);
% OptPrimalDual_counter   = zeros(3, Niter+1);
clear
load final-result.mat

%% counter calculation %%%%
comp_time_unit = 1;
comm_time_unit = 1;

% the first row is for total cost; the second row is for communication
% cost; the third row is for gradient computation cost.

% NEXT
NEXT_counter      = zeros(3, Niter+1);
NEXT_counter(2,:) = (0:Niter)*2*comm_time_unit;
NEXT_counter(3,:) = (0:Niter)*1*comp_time_unit;
NEXT_counter(1,:) = NEXT_counter(2,:) + NEXT_counter(3,:);

% EXTRA
EXTRA_counter          = zeros(3, Niter+1);
EXTRA_counter(2,2:end) = ((1:Niter)*2 - 1) * comm_time_unit;
EXTRA_counter(3,2:end) = (1:Niter)*1*comp_time_unit;
EXTRA_counter(1,:)     = EXTRA_counter(2,:) + EXTRA_counter(3,:);

% DPSGD
num_iter = ceil(2 * Niter);
DPSGD_counter      = zeros(3, num_iter+1);
DPSGD_counter(2,:) = (0:num_iter)*1*comm_time_unit;
DPSGD_counter(3,:) = (0:num_iter)*batch_size_portion*comp_time_unit;
DPSGD_counter(1,:) = DPSGD_counter(2,:) + DPSGD_counter(3,:);

% Acc_DNGD_NSC
Acc_DNGD_NSC_counter      = zeros(3, Niter+1); 
Acc_DNGD_NSC_counter(2,:) = (0:Niter)*3*comm_time_unit;
Acc_DNGD_NSC_counter(3,:) = (0:Niter)*1*comp_time_unit;
Acc_DNGD_NSC_counter(1,:) = Acc_DNGD_NSC_counter(2,:) + Acc_DNGD_NSC_counter(3,:);

% APM_C
L          = L_f;
Lap_eig    = sort(eig(Lap));
W          = eye(Num_Nodes) - Lap/Lap_eig(end);
W_eig      = sort(eig(W));

APM_C_counter           = zeros(3, Niter+1);
inner_consensus_factor  = zeros(1, Niter);
inner_consensus_const   = 0.2;

for k = 1:Niter
    inner_consensus_factor(k) = log(k) / sqrt(1-W_eig(end-1));
    T = ceil(inner_consensus_const * inner_consensus_factor(k));
    APM_C_counter(1,k+1) = APM_C_counter(1,k) + T*comm_time_unit + 1*comp_time_unit;
end
APM_C_counter(3,:) = (0:Niter)*1*comp_time_unit;
APM_C_counter(2,:) = APM_C_counter(1,:) - APM_C_counter(3,:);

% OptPrimalDual
K         = floor(Num_InnerConsensus/2);
OptPrimalDual_counter      = zeros(3, Niter+1);
OptPrimalDual_counter(2,:) = (0:Niter)*2*K*comm_time_unit;
OptPrimalDual_counter(3,:) = (0:Niter)*1*comp_time_unit;
OptPrimalDual_counter(1,:) = OptPrimalDual_counter(2,:) + OptPrimalDual_counter(3,:);

%% plot figures %%%%
figure
font_size = 16;
subplot(1,3,1)
fmean = semilogy(NEXT_counter(1,:), NEXT_err(:,1), EXTRA_counter(1,:), EXTRA_err(:,1),...
    Acc_DNGD_NSC_counter(1,:), Acc_DNGD_NSC_err(:,1),...
    APM_C_counter(1,:), APM_C_err(:,1), DPSGD_counter(1,:), DPSGD_err(:,1), ...
    OptPrimalDual_counter(1,:), OptPrimalDual_err(:,1));
hold on
set(fmean, 'linewidth', 2);
xlabel({'total cost'}, 'FontSize', font_size)
ylabel({'Bregman distance optimality gap'}, 'FontSize', font_size)
legend({'NEXT/DIGing', 'EXTRA', 'Acc-DNGD-NSC', 'APM-C', 'DPSGD', 'Our method'}, 'FontSize', font_size)
xlim([0 NEXT_counter(1,end)])
set(gca,'FontSize', font_size)

subplot(1,3,2)
fmean = semilogy(NEXT_counter(2,:), NEXT_err(:,1), EXTRA_counter(2,:), EXTRA_err(:,1),...
    Acc_DNGD_NSC_counter(2,:), Acc_DNGD_NSC_err(:,1),...
    APM_C_counter(2,:), APM_C_err(:,1), DPSGD_counter(2,:), DPSGD_err(:,1), ...
    OptPrimalDual_counter(2,:), OptPrimalDual_err(:,1));
hold on
set(fmean, 'linewidth', 2);
xlabel({'communication cost'}, 'FontSize', font_size)
% ylabel({'Bregman distance optimality gap'}, 'FontSize',20)
legend({'NEXT/DIGing', 'EXTRA', 'Acc-DNGD-NSC', 'APM-C', 'DPSGD', 'Our method'}, 'FontSize', font_size)
xlim([0 NEXT_counter(2,end)])
set(gca,'FontSize', font_size)

subplot(1,3,3)
fmean = semilogy(NEXT_counter(3,:), NEXT_err(:,1), EXTRA_counter(3,:), EXTRA_err(:,1),...
    Acc_DNGD_NSC_counter(3,:), Acc_DNGD_NSC_err(:,1),...
    APM_C_counter(3,:), APM_C_err(:,1), DPSGD_counter(3,:), DPSGD_err(:,1), ...
    OptPrimalDual_counter(3,:), OptPrimalDual_err(:,1));
hold on
% semilogy(OPTRA_counter(3,:), upperbound, 'ro');
set(fmean, 'linewidth', 2);
xlabel({'computation cost'}, 'FontSize', font_size)
% ylabel({'Bregman distance optimality gap'}, 'FontSize',20)
legend({'NEXT/DIGing', 'EXTRA', 'Acc-DNGD-NSC', 'APM-C', 'DPSGD', 'Our method'}, 'FontSize', font_size)
xlim([0 NEXT_counter(3,end)])
set(gca,'FontSize', font_size)


figure
subplot(1,3,1)
fmean = semilogy(NEXT_counter(1,:), NEXT_err(:,2), EXTRA_counter(1,:), EXTRA_err(:,2),...
    Acc_DNGD_NSC_counter(1,:), Acc_DNGD_NSC_err(:,2),...
    APM_C_counter(1,:), APM_C_err(:,2), DPSGD_counter(1,:), DPSGD_err(:,2), ...
    OptPrimalDual_counter(1,:), OptPrimalDual_err(:,2));
hold on
set(fmean, 'linewidth', 2);
xlabel({'total cost'}, 'FontSize', font_size)
ylabel({'function value optimality gap'}, 'FontSize', font_size)
legend({'NEXT/DIGing', 'EXTRA', 'Acc-DNGD-NSC', 'APM-C', 'DPSGD', 'Our method'}, 'FontSize', font_size)
xlim([0 NEXT_counter(1,end)])
set(gca,'FontSize', font_size)

subplot(1,3,2)
fmean = semilogy(NEXT_counter(2,:), NEXT_err(:,2), EXTRA_counter(2,:), EXTRA_err(:,2),...
    Acc_DNGD_NSC_counter(2,:), Acc_DNGD_NSC_err(:,2),...
    APM_C_counter(2,:), APM_C_err(:,2), DPSGD_counter(2,:), DPSGD_err(:,2), ...
    OptPrimalDual_counter(2,:), OptPrimalDual_err(:,2));
hold on
set(fmean, 'linewidth', 2);
xlabel({'communication cost'}, 'FontSize', font_size)
% ylabel({'function value optimality gap'}, 'FontSize',20)
legend({'NEXT/DIGing', 'EXTRA', 'Acc-DNGD-NSC', 'APM-C', 'DPSGD', 'Our method'}, 'FontSize', font_size)
xlim([0 NEXT_counter(2,end)])
set(gca,'FontSize', font_size)

subplot(1,3,3)
fmean = semilogy(NEXT_counter(3,:), NEXT_err(:,2), EXTRA_counter(3,:), EXTRA_err(:,2),...
    Acc_DNGD_NSC_counter(3,:), Acc_DNGD_NSC_err(:,2),...
    APM_C_counter(3,:), APM_C_err(:,2), DPSGD_counter(3,:), DPSGD_err(:,2), ...
    OptPrimalDual_counter(3,:), OptPrimalDual_err(:,2));
hold on
set(fmean, 'linewidth', 2);
xlabel({'computation cost'}, 'FontSize', font_size)
% ylabel({'function value optimality gap'}, 'FontSize',20)
legend({'NEXT/DIGing', 'EXTRA', 'Acc-DNGD-NSC', 'APM-C', 'DPSGD', 'Our method'}, 'FontSize', font_size)
xlim([0 NEXT_counter(3,end)])
set(gca,'FontSize', font_size)