function fm = fn_base2(ind)


e11 = cell(1,2);
e11{1,1} = 1;
e11{1,2} = [ind;1;1];

e12 = cell(1,2);
e12{1,1} = 1;
e12{1,2} = [ind;2;1];

e21 = cell(1,2);
e21{1,1} = 1;
e21{1,2} = [ind;1;2];

e22 = cell(1,2);
e22{1,1} = 1;
e22{1,2} = [ind;2;2];

fm = cell(2,2);
fm{1,1} = e11;
fm{1,2} = e12;
fm{2,1} = e21;
fm{2,2} = e22;
    
end