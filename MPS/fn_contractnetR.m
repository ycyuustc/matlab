function [R]=fn_contractnetR(mps,m)

R=1;
N=length(mps);

for k=N:-1:(m+1)
   
    R=fn_contract(conj(mps{k}),3,2,R,2,1);
    R=fn_contract(mps{k},3,[2,3],R,3,[3,2]);
    R=permute(R,[2,1]);
    
end

end