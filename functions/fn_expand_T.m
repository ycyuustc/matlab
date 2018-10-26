function [fc, fc1] = fn_expand_T(T, ind1, ind2)
    
fc = cell(16,2);
point = 1;
for alpha1=1:2
    for alpha2=1:2
        for beta1=1:2
            for beta2=1:2
                
                para = T(alpha1,alpha2,beta1,beta2);
                matrix = [ind1,ind2;beta1,beta2;alpha1,alpha2];
                fc{point,1} = para;
                fc{point,2} = matrix;
                point = point + 1;
                
            end
        end
    end
end

fc1 = fc;
fc = fn_fermion_comb(fc);
    
end