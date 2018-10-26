function fc = fn_fermion_comb(fc_input)
 
fc = fc_input;
point = 1;

while 1
    [num_base, ~] = size(fc);
    if point>=num_base
        break;
    end
    
    one_hot = zeros(1,num_base);
    select_base = fc{point,2};
    
    sum_para = fc{point,1};
    for i=(point+1):num_base
        if norm(size(select_base) - size(fc{i,2}))<1e-1
            if norm(select_base - fc{i,2})<1e-2
                one_hot(i) = 1;
                sum_para = sum_para + fc{i,1};
            end
        end
    end
    
    if isa(sum_para,'sym')
        fc{point,1} = simplify(sum_para);
    else
        fc{point,1} = sum_para;
    end
    
    fc(one_hot==1,:) = [];
    point = point + 1;
    
end

[num_base, ~] = size(fc);
v_onehot = zeros(1,num_base);
for i=1:num_base  
    if isa(fc{i,1},'double')
        if abs(fc{i,1})<1e-8
            v_onehot(i) = 1;
        end 
    end
    if isa(fc{i,1},'sym')
        if fc{i,1} == '0'
            v_onehot(i) = 1;
        end
    end
end

fc(v_onehot==1,:) = [];

end