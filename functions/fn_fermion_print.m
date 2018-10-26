function fc_string = fn_fermion_print(fc)

[num2,~] = size(fc);
if num2 == 0
    fc_string ='0';
    return;
end
cpara = cell(1,num2);
cbase = cell(1,num2);
for k=1:num2
    
    matrix = fc{k,2};
    [~,num_b] = size(matrix);
    temp = [];
    for k2 = 1:num_b
        
        index = matrix(1,k2);
        ind1 = matrix(2,k2);
        ind2 = matrix(3,k2);
              
        add_on = ['[','e',num2str(index),'_',num2str(ind2),num2str(ind1),']'];
        temp = strcat(temp,add_on);
        %            temp = strcat(temp,'e');
        %            temp = strcat(temp,num2str(k2));
        %            temp = strcat(temp,'_');
        %            temp = strcat(temp,num2str(ind2));
        %            temp = strcat(temp,num2str(ind1));
        
    end
    
    cbase{1,k} = temp;
    para = fc{k,1};
    if isa(para,'double')
        para = num2str(para);
    else
        para = char(para);
    end
    cpara{1,k} = para;
    
    
end

fc_string = [];
for k=1:num2
    
    para = cpara{1,k};
    tbase = cbase{1,k};
    
    if para(1)~='-'
        para = strcat('+',para);
    end
    %           if para(2)=='1' && length(para)>2
    %               para(2) = [];
    %           end
    
    fc_string = strcat(fc_string,para);
    fc_string = strcat(fc_string,tbase);
    fc_string = strcat(fc_string,'\n');
    
end


if fc_string(1)=='+'
    fc_string(1) = [];
end



end