function fct = fn_fermion_convert(fc, ind)
    
   [num2,~] = size(fc);
   add_on = [];
   for k=1:num2
      
       matrix = fc{k,2};
       tv = matrix(1,:);
       index = find(tv==ind);
       index = index(1);
       tv2 = matrix(:,index);
       if tv2(2) == 2 && tv2(3) ==2
          
           tv2(2) = 1;
           tv2(3) = 1;
           matrix(:,index) = tv2;
           fc{k,2} = matrix;
           fc{k,1} = -fc{k,1};
           
           temp = cell(1,2);
           temp{1,1} = -fc{k,1};
           matrix(:,index) = [];
           temp{1,2} = matrix;
           add_on = fn_fermion_add(add_on,temp);
           
       end
       
   end
   
   fct = fn_fermion_add(fc, add_on);
   fct = fn_fermion_comb(fct);
    
end
