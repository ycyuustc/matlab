function [result,re_sign] = fn_fermion_base_product(fm1, fm2)

if isempty(fm1)
    result = fm2;
    re_sign = 1;
    return;
end

if isempty(fm2)
    result = fm1;
    re_sign = 1;
    return;
end

tv1 = fm1(1,:)*10 + 1;
tv2 = fm2(1,:)*10 + 2;
fm1(1,:) = tv1;
fm2(1,:) = tv2;

fm = [fm1, fm2];
tv = fm(1,:);
num_base = length(tv);
sign_sum = 0;
for i1=1:num_base
    for i2=(i1+1):num_base
        if fm(1,i1)>fm(1,i2)
            sign_sum = sign_sum + (fm(2,i1)+fm(3,i1))*(fm(2,i2)+fm(3,i2)); 
        end
    end
end
re_sign = 1 - 2*mod(sign_sum,2);

[~,sort_tv] = sort(tv);
fm = fm(:,sort_tv);

result = zeros(size(fm));
point = 1;
point_w = 1;
while 1
    
    p1 = point;
    p2 = point + 1;
    ind1 = fm(1,p1);
    tv1 = fm(2:3,p1);
    if p2<=num_base
        ind2 = fm(1,p2);
        tv2 = fm(2:3,p2);
    else
        ind2 = ind1 + 5;
    end
    
    if ind2 - ind1 <2
        if tv1(1) ~= tv2(2)
            result = [];
            re_sign = 0;
            break;
        end      
        tv = [tv2(1); tv1(2)];
        result(:,point_w) = [floor(ind1/10.0); tv];
        point_w = point_w + 1;
        point = point + 2;
    else
        result(:,point_w) = fm(:, point);
        result(1,point_w) = floor(result(1, point_w)/10);
        point = point + 1;
        point_w = point_w + 1;
    end
    
    if point>num_base
        result = result(:,1:(point_w-1));
        break;
    end
    
end

    
end