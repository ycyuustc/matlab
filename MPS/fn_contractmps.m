function result=fn_contractmps(mps)
% contract a simple mps

N=length(mps);

L=1;

for k=1:N
    
    L=fn_contract(L,2,1,conj(mps{k}),3,1);
    L=fn_contract(L,3,[1,3],mps{k},3,[1,3]);

end

result=L;