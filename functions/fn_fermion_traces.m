function fct = fn_fermion_trace(fc, ind)

fct = [];
for p=1:2
    for q=1:2
        
        pleft = cell(1,2);
        pright = cell(1,2);
        pleft{1,1} = 1;
        pright{1,1} = 1;
        pleft{1,2} = [ind;q;p];
        pright{1,2} = [ind;p;q];
        temp = fn_fermion_product(pleft, fc);
        temp = fn_fermion_product(temp,pright);
        fct = fn_fermion_add(fct,temp);
        
    end
end

fct = fn_fermion_comb(fct);
    
end