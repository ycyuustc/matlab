function [B,K]=fn_reduceD2_onesite(A,X,Cleft,Cright)

Cleft=fn_contract(Cleft,3,3,A,3,1);
Cleft=fn_contract(Cleft,4,[2,4],X,4,[1,4]);
B=fn_contract(Cleft,4,[3,2],Cright,3,[2,3]);
B=permute(B,[1,3,2]);
b=reshape(B,[numel(B),1]);
K=-b'*b;

end