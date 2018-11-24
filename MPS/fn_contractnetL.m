function [L]=fn_contractnetL(mps,m)

L=1;
% N=length(mps);

for k=1:(m-1)
   
    L=fn_contract(L,2,1,conj(mps{k}),3,1);
    L=fn_contract(L,3,[1,3],mps{k},3,[1,3]);
    
end

end