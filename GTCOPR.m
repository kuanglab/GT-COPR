function A = GTCOPR(W,R,opts)

% initialize the CPD form
net_num = length(W);
nvec = zeros(1,net_num);
A = cell(1,net_num);
for netid = 1:length(W)
rng(1);
nvec(netid) = size(W{netid},1);
A{netid}=rand(nvec(netid),opts.rank_k);
end

% compute the auxiliary variables
[R,fac] = tensor_R(R,net_num,nvec,opts.normR);
[X,Y,ZW,ZD,MW,MD,W,D] = buil_AUX(A,W,net_num,opts.normW);

Uori = cell(1,net_num);
for i=1:net_num
    Uori{i} = sum(A{i},1);
end

% multiplicative updating
for iter=1:opts.maxIter
    Aold = A;
    for netid=1:net_num
        [num,denom] = update_collap(A,R,X,Y,nvec,Uori,fac,net_num,netid);     
        [num_graph,denom_graph] = update_graphReg(A,net_num,netid,X,ZW,ZD,MW,MD,opts.graphRegType);
        num=num+opts.lambda*num_graph+10^-10;
        denom=denom+opts.lambda*denom_graph+opts.beta*A{netid}+10^-10;
        A{netid}=A{netid}.*(num./denom).^(opts.stepsize);
        [X,Y,ZW,ZD,MW,MD] = update_AUX(X,Y,ZW,ZD,MW,MD,A,W,D,netid);
        Uori{netid} = sum(A{netid},1);
    end
res = compute_res(A,Aold);    
if mod(iter,1)==0
   disp(['res: ',num2str(res)]);
end
if res<opts.thre
    break;
end

end
end

function res = compute_res(A,Aold)
res_num = 0;
res_denom = 0;
for i = 1:length(A)
    res_num = res_num + sum(sum((A{i}-Aold{i}).^2));
    res_denom = res_denom + sum(sum(Aold{i}.^2));
end
res=sqrt(res_num/res_denom);
end
