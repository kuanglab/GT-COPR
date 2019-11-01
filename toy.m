clear all;clc
opts.lambda = 10^-2; % hyperparameter of graph regularization
opts.beta = 10^-3; % hyperparameter of L-2 regularization
opts.thre = 0.01;  % stopping condition
opts.stepsize = 1; % step size of multiplicative updating
opts.normR = 1; % normalize the pairwise relations? (1: Yes, 0: No)
opts.normW = 0; % normalize the knowledge graphs? (1: Yes, 0: No)
opts.rank_k = 10; % CPD rank
opts.maxIter = 300; % the maximum # of iterations
opts.graphRegType = 'SPG'; % the type of product graph ('SPG': strong, 'TPG': tensor, 'CPG': Cartesian)
nvec=[100,90,70,50]; % sizes of graphs
net_num=length(nvec); % number of graphs
W=cell(net_num,1); % store the graph adjacency matrices in a net_num x 1 cell array
for i=1:net_num
    rng(1);
    W{i} = rand(nvec(i),nvec(i));
    W{i} = W{i} + W{i}';
    W{i} = W{i} - diag(diag(W{i}));
end
R = cell(net_num,net_num);
for i=1:net_num-1
    for j=i+1:net_num
        R{i,j} = rand(nvec(i),nvec(j)); % store the pairwise relations in the upper triangular of a net_num x net_num cell array
    end
end
Ares = GTCOPR(W,R,opts); % run GT-COPR algorithm and save the learned factor matrices in a net_num x 1 cell array
