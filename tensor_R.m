function [R,fac] = tensor_R(R,net_num,nvec,normR)
if  normR==1
    R = norm_R(R,net_num);
end
Rt = cell(net_num,net_num);
fac = Rt;
dims_ori = ones(size(nvec));
fullset = 1:net_num;
for i=1:net_num-1
    for j=i+1:net_num
        [ia,ib,ic] = find(R{i,j});
        subs = ones(length(ia),net_num);
        subs(:,[i,j]) = [ia,ib];
        dims = dims_ori;
        dims(i) = nvec(i);
        dims(j) = nvec(j);
        Rt{i,j} = sptensor(subs,ic,dims);
        fac{i,j} = prod(nvec(setdiff(fullset,[i,j])));
        Rt{i,j} = Rt{i,j}*fac{i,j};
    end
end
R = Rt;
end

function R = norm_R(R,net_num)
 for i=1:net_num-1
     for j=i+1:net_num
         mat1=sqrt(repmat(sum(R{i,j},1),size(R{i,j},1),1));
         mat2=sqrt(repmat(sum(R{i,j},2),1,size(R{i,j},2)));
         R{i,j}(mat1~=0)=R{i,j}(mat1~=0)./mat1(mat1~=0);
         R{i,j}(mat2~=0)=R{i,j}(mat2~=0)./mat2(mat2~=0);
         R{i,j} = (R{i,j} - min(R{i,j}(:)))/(max(R{i,j}(:))-min(R{i,j}(:)));
     end
 end
end