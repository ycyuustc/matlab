function [result, matrix] = fn_cell_kron(cA,cB)
   
    [a1,a2]=size(cA);
    [b1,b2]=size(cB);
    
    result = cell(a1,b1,a2,b2);
    for ia1=1:a1
        for ib1=1:b1
            for ia2=1:a2
                for ib2=1:b2
                    result{ia1,ib1,ia2,ib2} = simplify(cA{ia1,ia2}*cB{ib1,ib2});
                end
            end
        end
    end
    
    matrix = reshape(permute(result,[2,1,4,3]),[4,4]);
    
end