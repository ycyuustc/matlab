function [mpsA,mpsB,cs]=fn_prepare_mps_ns(mpsA,mpsB)

N=length(mpsA);
cs = cell(1,N);
cs{N} = 1;

for i=N:-1:2
    
    [mpsA{i},mpsB{i},A_out,B_out,AB_s]...
        =fn_prepare_onesite_rl_ns(mpsA{i},mpsB{i},cs{i});
    tcell=fn_contract(mpsA{i-1},3,2,A_out,2,1);
    mpsA{i-1}=permute(tcell,[1,3,2]);
    tcell=fn_contract(mpsB{i-1},3,2,B_out,2,1);
    mpsB{i-1}=permute(tcell,[1,3,2]);
    
    cs{i-1} = AB_s;
    
end

end