function [num,denom] = update_collap(A,R,X,Y,nvec,Uori,fac,net_num,netid)
ker = cell(net_num,net_num);
num = zeros(size(A{netid}));
rank_k = size(A{1},2);
denom1 = zeros(rank_k,rank_k);
denom2 = denom1;
for i=1:net_num-1
    for j=i+1:net_num
        U = Uori;
        U{i} = A{i};
        U{j} = A{j};
        ker{i,j} = mttkrp(R{i,j},U,netid)/(fac{i,j}^2);
        if (i ~= netid) && (j ~= netid)
            ker{i,j} = repmat(ker{i,j},nvec(netid),1);
        end

        mat = ones(rank_k,rank_k);
        for id =1:net_num
            if id ~= netid
                mat_tmp = X{id};
                if (i ~= id) && (j ~= id)
                    mat_tmp = Y{id};
                end    
                mat = mat.*mat_tmp;
            end
        end
        mat = mat/(fac{i,j}^2);
        if (i ~= netid) && (j ~= netid)
            denom2 = denom2 + mat;

        else
            denom1 = denom1 + mat;
        end

        num = num + ker{i,j};  
    end
end
denom = A{netid}*denom1 + repmat(sum(A{netid},1)*denom2,nvec(netid),1);
end