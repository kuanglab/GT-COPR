function [num_graph,denom_graph] = update_graphReg(A,net_num,netid,X,ZW,ZD,MW,MD,graphRegType)
if strcmp(graphRegType,'CPG')==1
    [num_graph,denom_graph] = CPG(A,net_num,netid,X,ZW,ZD,MW,MD);
elseif strcmp(graphRegType,'TPG')==1
    [num_graph,denom_graph] = TPG(A,net_num,netid,ZW,ZD,MW,MD);
else
    [num_graph,denom_graph] = SPG(A,net_num,netid,X,ZW,ZD,MW,MD);
end

end

function [num_CPG,denom_CPG] = CPG(A,net_num,netid,X,ZW,ZD,MW,MD)
rank_k = size(A{1},2);
num_G = zeros(rank_k,rank_k);
denom_G = num_G;
numdenom = num_G;
for i=1:net_num-1
    for j=i+1:net_num                
        mat_num = ones(rank_k,rank_k);
        mat_denom = ones(rank_k,rank_k);
        for id =1:net_num
            if id ~= netid
                mat_tmpnum = X{id};
                mat_tmpdenom = X{id};
                if (i ~= id) && (j ~= id)
                    mat_tmpnum = MW{id};
                    mat_tmpdenom = MD{id};
                end    
                mat_num = mat_num.*mat_tmpnum;
                mat_denom = mat_denom.*mat_tmpdenom;
            end
        end
        if (i ~= netid) && (j ~= netid)
            numdenom = numdenom + mat_num;
        else
            num_G = num_G + mat_num;
            denom_G = denom_G + mat_denom;
        end

    end
end
num_CPG=A{netid}*num_G+ZW{netid}*numdenom;
denom_CPG=A{netid}*denom_G+ZD{netid}*numdenom;
end


function [num_TPG,denom_TPG] = TPG(A,net_num,netid,ZW,ZD,MW,MD)
rank_k = size(A{1},2);
mat_num = ones(rank_k,rank_k);
mat_denom = ones(rank_k,rank_k);
for id =1:net_num
    if id ~= netid
        mat_tmpnum = MW{id};
        mat_tmpdenom = MD{id};  
        mat_num = mat_num.*mat_tmpnum;
        mat_denom = mat_denom.*mat_tmpdenom;
    end
end
num_TPG = ZW{netid}*mat_num;
denom_TPG = ZD{netid}*mat_denom;
end


function [num_SPG,denom_SPG] = SPG(A,net_num,netid,X,ZW,ZD,MW,MD)
[num_CPG,denom_CPG] = CPG(A,net_num,netid,X,ZW,ZD,MW,MD);
[num_TPG,denom_TPG] = TPG(A,net_num,netid,ZW,ZD,MW,MD);
num_SPG = num_CPG + num_TPG;
denom_SPG = denom_CPG + denom_TPG;
end