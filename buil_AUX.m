function [X,Y,ZW,ZD,MW,MD,W,D] = buil_AUX(A,W,net_num,normW)
X = cell(net_num,1);
Y = X;
ZW = X;
ZD = X;
MW = X;
MD = X;
if normW == 1
    for netid=1:net_num 
        d=sum(W{netid},2);
        d(d~=0)=(d(d~=0)).^-(0.5);
        W{netid}=sparse(diag(d))*sparse(W{netid})*sparse(diag(d));
    end
end
D=cell(net_num,1);
for netid=1:net_num
    D{netid}=diag(sum(W{netid},2));
end 

for i=1:net_num
    [X,Y,ZW,ZD,MW,MD] = update_AUX(X,Y,ZW,ZD,MW,MD,A,W,D,i);
end
end

