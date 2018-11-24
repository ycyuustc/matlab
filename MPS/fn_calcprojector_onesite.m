function [P]=fn_calcprojector_onesite(B,Cleft,Cright)

y=fn_contract(Cleft,3,3,B,3,1);
y=fn_contract(y,4,[2,3],Cright,3,[2,3]);
y=permute(y,[1,3,2]);
y=reshape(y,[numel(y),1]);
Q=orth([y,eye(size(y,1))]);

P=Q(:,2:end);

end