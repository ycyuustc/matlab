function fc = fn_to_vacuum(fc)
   
    [num2,~] = size(fc);
    for k = 1:num2
       mat_base = fc{k,2};
       [~, num_base] = size(mat_base);
       flag = 1;
       for k2 = 1:num_base
          if mat_base(2,k2) == 2
              flag = 0;
          end
       end
       fc{k,1} = flag*fc{k,1};
    end
    fc = fn_fermion_comb(fc);
    
end