function fc = fn_fermion_minus(fc1,fc2)
   
    [num2,~] = size(fc2);
    for k=1:num2
        fc2{k,1} = - fc2{k,1};
    end
    
    fc = [fc1;fc2];
    fc = fn_fermion_comb(fc);
    
end