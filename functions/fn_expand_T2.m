function fm = fn_expand_T2(T,ind1)
    
    fm = cell(2,2);
    
    for alpha1=1:2
        for beta1=1:2
            
            fc = cell(4,2);
            point = 1;
            for alpha2=1:2
                for beta2=1:2
                    
                    para = T(alpha1,alpha2,beta1,beta2);
                    fc{point,1} = para;
                    fc{point,2} = [ind1;beta2;alpha2];
                    point = point + 1;
                    
                end
            end
            
            fc = fn_fermion_comb(fc);
            fm{alpha1,beta1} = fc;
            
        end
    end
    
end