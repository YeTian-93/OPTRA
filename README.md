# OPTRA
  
This is the MATLAB code for simulating OPTRA.  The simulation consists of two problems--decentralized linear regression and decentralized logistic regression.  We compare the proposed rate-optimal algorithm–OPTRA–with existing accelerated ones designed for convex smooth problems, namely: Acc-DNGD-NSC (Qu and Li, 2017a) and APM-C (Li et al., 2018). We also included non-accelerated schemes that perform quite well in practice, namely: i) the gradient tracking method, NEXT/DIGing (Di Lorenzo and Scutari, 2016; Nedich et al., 2017); ii) the primal-dual method, EXTRA (Shi et al., 2015); and iii) the decentralized stochastic gradient method, DPSGD (Lian et al., 2017).  The detailed setting of the simulation can be found in our paper [1].

Each folder contains codes for the specific problem.  Inside each folder, the main function is main.m.  Please contact me via _tian110@purdue.edu_ for any question.  Enjoy!

## Reference:
[1] Jinming Xu, Ye Tian, Ying Sun, and Gesualdo Scutari. "Accelerated Primal-Dual Algorithms for Distributed Smooth Convex Optimization over Networks." In The 23rd International Conference on Artificial Intelligence and Statistics (AISTATS), 2020.
