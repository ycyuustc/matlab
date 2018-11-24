function [Cright]=fn_updateCright(Cright,B,X,A)

if isempty(X)
   d=size(B,3);
   X=reshape(eye(d),[1,1,d,d]);
end

Cright=fn_contract(A,3,2,Cright,3,3);
Cright=fn_contract(X,4,[2,4],Cright,4,[4,2]);
Cright=fn_contract(conj(B),3,[2,3],Cright,4,[4,2]);

end