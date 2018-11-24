function result=fn_contractmps2(mps1,mps2)
% contract a simple mps

N=length(mps1);

L=1;

for k=1:N
    
    L=fn_contract(L,2,1,conj(mps2{k}),3,1);
    L=fn_contract(L,3,[1,3],mps1{k},3,[1,3]);

end

result=L;