1) main.m is the main function;  
2) EXTRA.m, NEXT.m, Acc_DNGD_NSC.m, APM_C.m, OPTRA.m are the five algorithms we simulate;
3) LeastSquare_Data_Gen.m is to generated the problem data;
4) Network_Laplacian.m is to generate a Laplacian matrix of a connected graph;
5) AccGossip.m is to implement the accelerated gossip by Chebyshev acceleration;
6) data.mat and laplacian.mat provide instances of the problem data and the Laplacian matrix; 
7) evaluate.m applies the Bregman distance optimality measure and the function value optimality measure on an input iteration;
8) result.mat is the result data obtained by running main.m 

The user can run main.m directly to simulate and generate the result.