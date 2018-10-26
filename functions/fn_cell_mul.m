function result = fn_cell_mul(cA,cB)
   
    [a1,a2]=size(cA);
    [b1,b2]=size(cB);
    
    if a2~=b1
        error('dimensions are not feat');
    end
    
    result=cell(a1,b2);
    
    
    for x=1:a1
        for y=1:b2
            temp = 0;
            for k = 1:a2
                temp = temp + cA{x,k}*cB{k,y};
            end
            result{x,y}=simplify(temp);
        end
    end
       
end