function [Cleft]=fn_updateCleft(Cleft,B,X,A)

if isempty(X)
   d=size(B,3);
   X=reshape(eye(d),[1,1,d,d]);
end

Cleft=fn_contract(A,3,1,Cleft,3,3);
Cleft=fn_contract(X,4,[1,4],Cleft,4,[4,2]);
Cleft=fn_contract(conj(B),3,[1,3],Cleft,4,[4,2]);

end

