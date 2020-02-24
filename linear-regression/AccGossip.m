function newLap = AccGossip(Lap, K, c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code implement the accelerated gossip by Chebyshev acceleration 
% first introduced in [1] (Algorithm 2).  Our implementation produces an 
% equivalent new Laplacian matrix after the accelerated gossip.
% ---------input----------
% Lap:     a Laplacian matrix;
% K:       number of steps for inner consensus;
% c:       constant c_1 specified in Algorithm 2 in our paper
% ---------output---------
% newLap:  the equivalent new Laplacian matrix
% --------reference-------
% [1] Seaman, Kevin, Francis Bach, Sébastien Bubeck, Yin Tat Lee, and 
%     Laurent Massoulié. "Optimal algorithms for smooth and strongly convex 
%     distributed optimization in networks." In Proceedings of the 34th 
%     International Conference on Machine Learning-Volume 70, pp. 3027-3036. 
%     JMLR. org, 2017.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Num_Nodes
a = 1;   new_a = c;
Z = eye(Num_Nodes);   new_Z = c*(eye(Num_Nodes) - Lap);
for k = 1:K-1
    temp_a = 2*c*new_a - a;
    temp_Z = 2*c*(eye(Num_Nodes) - Lap)*new_Z - Z;
    a = new_a; Z = new_Z;
    new_a = temp_a; new_Z = temp_Z;
end
newLap = eye(Num_Nodes) - new_Z/new_a;