function result = fn_superize(ca)
    
   [a1,a2] = size(ca);
   
   if a1~=4 || a2~=4
       error('dimension is not equal to 4');
   end
   
   result = ca;
   result{2,3} = -result{2,3};
   result{2,4} = -result{2,4};
   result{4,1} = -result{4,1};
   result{4,2} = -result{4,2};
    
end