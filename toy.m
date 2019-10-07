clear all;clc
opts.lambda = 10^-2;
opts.beta = 10^-3;
opts.thre = 0.01;
opts.stepsize = 1;
opts.normR = 1;
opts.normW = 0;
opts.rank_k = 10;
opts.maxIter = 300;
opts.graphRegType = 'SPG';
nvec=[100,90,70,50];
net_num=length(nvec);
W=cell(net_num,1);
for i=1:net_num
    rng(1);
    W{i} = rand(nvec(i),nvec(i));
    W{i} = W{i} + W{i}';
    W{i} = W{i} - diag(diag(W{i}));
end
R = cell(net_num,net_num);
for i=1:net_num-1
    for j=i+1:net_num
        R{i,j} = rand(nvec(i),nvec(j)); 
    end
end
Ares = GTCOPR(W,R,opts);
