function res = mcon(cTensor,cIndex)

numTen = length(cIndex);
numIndex = 0;
numIndexOut = 0;
for i=1:numTen
    
    numIndex = numIndex + length(cIndex{i});
    numIndexOut = numIndexOut + sum(cIndex{i}<0);
    %     disp(numIndexOut);
    %     disp(numIndex);
    
end

numIndexIn = (numIndex - numIndexOut)/2;

v_order(netcon(cIndex,0,2,1,1)) = 1:numIndexIn;
c_order = cell(1,numTen);
for i=1:numTen
    
    tv = cIndex{i};
    for j = 1:length(tv)
        if tv(j)>0
            tv(j) = v_order(tv(j));
        end
    end
    
    c_order{i} = tv;
    
end

res = ncon(cTensor,c_order);

end