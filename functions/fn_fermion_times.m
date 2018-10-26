function fc = fn_fermion_times(lambda, fc2)
   
    [num2,~] = size(fc2);
    for k=1:num2
        fc2{k,1} = lambda * fc2{k,1};
    end
    
    fc = fn_fermion_comb(fc2);
    
end