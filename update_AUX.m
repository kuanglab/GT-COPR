function [X,Y,ZW,ZD,MW,MD] = update_AUX(X,Y,ZW,ZD,MW,MD,A,W,D,i)
X{i} = A{i}'*A{i};
Y{i} = sum(A{i},1)'*sum(A{i},1);
ZW{i} = W{i}*A{i};
ZD{i} = diag(D{i}).*A{i};
MW{i} = A{i}'*ZW{i};
MD{i} = A{i}'*ZD{i}; 
end