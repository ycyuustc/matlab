function [fc,fc_un] = fn_fermion_product(fc1, fc2)
    
    fc1 = fn_fermion_comb(fc1);
    fc2 = fn_fermion_comb(fc2);
   
    [num1,~] = size(fc1);
    [num2,~] = size(fc2);
    
    fc = cell(num1*num2,2);
    
    v_onehot = zeros(1,num1*num2);
    point = 1;
    for i1=1:num1
        for i2=1:num2
           
            para1 = fc1{i1,1};
            base1 = fc1{i1,2};
            para2 = fc2{i2,1};
            base2 = fc2{i2,2};
            
            [base12, sign12] = fn_fermion_base_product(base1, base2);
            temp = para1*para2*sign12;
            if isa(temp,'sym')
                para12 = simplify(temp);
            else
                para12 = temp;
            end
                   
            if para12 == 0
                v_onehot(point) = 1;
            end
            
            fc{point,1} = para12;
            fc{point,2} = base12;
            point = point + 1;
            
        end
    end
    
    fc_un = fc;
    fc(v_onehot==1,:) = [];
       
    fc = fn_fermion_comb(fc);
    
end