function fm = fn_expand_mul(T1, T2)
    
    fm = cell(2,2);
    for ind1=1:2
        for ind2=1:2
            
            tc1 = fn_fermion_product(T1{ind1,1},T2{1,ind2});
            tc2 = fn_fermion_product(T1{ind1,2},T2{2,ind2});
            
            tc = [tc1;tc2];
            fm{ind1,ind2} = fn_fermion_comb(tc);
            
        end
    end
    
    
    
end