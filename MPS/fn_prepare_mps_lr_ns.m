function [mpsA,mpsB,cs]=fn_prepare_mps_lr_ns(mpsA,mpsB)

N=length(mpsA);
cs = cell(1,N);
cs{1} = 1;

for i=1:(N-1)
    
    [mpsA{i},mpsB{i},A_out,B_out,AB_s]...
        =fn_prepare_onesite_lr_ns(mpsA{i},mpsB{i},cs{i});
    mpsA{i+1}=fn_contract(A_out,2,2,mpsA{i+1},3,1);

    mpsB{i+1}=fn_contract(B_out,2,2,mpsB{i+1},3,1); 
    
    cs{i+1} = AB_s;
    
end

end