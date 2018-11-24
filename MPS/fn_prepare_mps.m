function [mps]=fn_prepare_mps(mps)

N=length(mps);

for i=N:-1:2
   
    [mps{i},U]=fn_prepare_onesite_rl(mps{i});
    tcell=fn_contract(mps{i-1},3,2,U,2,1);
    mps{i-1}=permute(tcell,[1,3,2]);
    
end

end