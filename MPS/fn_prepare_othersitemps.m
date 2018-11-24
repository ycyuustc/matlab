function mps=fn_prepare_othersitemps(mps,m)

N=length(mps);

for k=N:-1:(m+1)
   
    [mps{k},U]=fn_prepare_onesite_rl(mps{k});
    mps{k-1}=fn_contract(mps{k-1},3,2,U,2,1);
    mps{k-1}=permute(mps{k-1},[1,3,2]);
    
end

for k=1:(m-1)
   
    [mps{k},U]=fn_prepare_onesite_lr(mps{k});
    mps{k+1}=fn_contract(U,2,2,mps{k+1},3,1);
    
end



end