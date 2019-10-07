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

% function [num,denom] = update_collap(A,R,X,Y,nvec,Uori,fac,net_num,netid)
% ker = cell(net_num,net_num);
% num = zeros(size(A{netid}));
% denom1 = num;
% denom2 = num;
% for i=1:net_num-1
%     for j=i+1:net_num
%         U = Uori;
%         U{i} = A{i};
%         U{j} = A{j};
%         ker{i,j} = mttkrp(R{i,j},U,netid)/(fac{i,j}^2);
%         if (i ~= netid) && (j ~= netid)
%             ker{i,j} = ones(nvec(netid),1)*ker{i,j};
%         end
% 
%         mat = ones(rank_k,rank_k);
%         for id =1:net_num
%             if id ~= netid
%                 mat_tmp = X{id};
%                 if (i ~= id) && (j ~= id)
%                     mat_tmp = Y{id};
%                 end    
%                 mat = mat.*mat_tmp;
%             end
%         end
%         mat = mat/(fac{i,j}^2);
%         if (i ~= netid) && (j ~= netid)
%             denom2 = denom2 + mat;
%         else
%             denom1 = denom1 + mat;
%         end
% 
%         num = num + ker{i,j};
%         disp('kronsum regularization');    
%     end
% end
% denom = A{netid}*denom1 + repmat(sum(A{netid},1)*denom2,nvec(netid),1);
% end



% function [num_kronsum,denom_kronsum] = update_graphReg(A,net_num,netid,rank_k,X,ZW,ZD,MW,MD)
% num_G = zeros(size(A{netid}));
% denom_G = num_G;
% numdenom = num_G;
% for i=1:net_num-1
%     for j=i+1:net_num                
%         mat_num = ones(rank_k,rank_k);
%         mat_denom = ones(rank_k,rank_k);
%         for id =1:net_num
%             if id ~= netid
%                 mat_tmpnum = X{id};
%                 mat_tmpdenom = X{id};
%                 if (i ~= id) && (j ~= id)
%                     mat_tmpnum = MW{id};
%                     mat_tmpdenom = MD{id};
%                 end    
%                 mat_num = mat_num.*mat_tmpnum;
%                 mat_denom = mat_denom.*mat_tmpdenom;
%             end
%         end
%         if (i ~= netid) && (j ~= netid)
%             numdenom = numdenom + mat_num;
%         else
%             num_G = num_G + mat_num;
%             denom_G = denom_G + mat_denom;
%         end
% 
%     end
% end
% num_kronsum=A{netid}*num_G+ZW{netid}*numdenom;
% denom_kronsum=A{netid}*denom_G+ZD{netid}*numdenom;
% end


% function [R,fac] = tensor_R(R,net_num,normR)
% if  normR==1
%     R = norm_R(R,net_num);
% end
% Rt = cell(net_num,net_num);
% fac = Rt;
% for i=1:net_num-1
%     for j=i+1:net_num
%         [ia,ib,ic] = find(R{i,j});
%         subs = ones(length(ia),net_num);
%         subs(:,[i,j]) = [ia,ib];
%         Rt{i,j} = sptensor(subs,ic,nvec);
%         fac{i,j} = prod(setdiff([1:net_num],[i,j]));
%         Rt{i,j} = Rt{i,j}*fac{i,j};
%     end
% end
% R = Rt;
% end
% 
% function R = norm_R(R,net_num)
%  for i=1:net_num-1
%      for j=i+1:net_num
%          mat1=sqrt(repmat(sum(R{i,j},1),size(R{i,j},1),1));
%          mat2=sqrt(repmat(sum(R{i,j},2),1,size(R{i,j},2)));
%          R{i,j}(mat1~=0)=R{i,j}(mat1~=0)./mat1(mat1~=0);
%          R{i,j}(mat2~=0)=R{i,j}(mat2~=0)./mat2(mat2~=0);
%          R{i,j} = (R{i,j} - min(R{i,j}(:)))/(max(R{i,j}(:))-min(R{i,j}(:)));
%      end
%  end
% end
% function [X,Y,ZW,ZD,MW,MD,W,D] = buil_AUX(A,W,net_num,normW)
% X = cell(net_num,1);
% Y = X;
% ZW = X;
% ZD = X;
% MW = X;
% MD = X;
% if normW == 1
%     for netid=1:net_num 
%         d=sum(W{netid},2);
%         d(d~=0)=(d(d~=0)).^-(0.5);
%         W{netid}=sparse(diag(d))*sparse(W{netid})*sparse(diag(d));
%     end
%    
% end
% D=cell(net_num,1);
% for netid=1:net_num
%     D{netid}=diag(sum(W{netid},2));
% end 
% 
% for i=1:net_num
%     [X,Y,ZW,ZD,MW,MD] = update_AUX(X,Y,ZW,ZD,MW,MD,A,W,D,i);
% end
% end
% 
% function [X,Y,ZW,ZD,MW,MD] = update_AUX(X,Y,ZW,ZD,MW,MD,A,W,D,i)
% X{i} = A{i}'*A{i};
% Y{i} = sum(A{i},1)'*sum(A{i},1);
% ZW{i} = W{i}*A{i};
% ZD{i} = diag(D{i}).*A{i};
% MW{i} = A{i}'*ZW{i};
% MD{i} = A{i}'*ZD{i}; 
% end
