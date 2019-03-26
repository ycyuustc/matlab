num_v=2312;
Num_w=1;
g=0;
phi=0.3;
beth=(sqrt(5.0)-1.0)/2;
v_t = cos(2.0*pi*beth*(1:1:30)+phi);

w =linspace(3,3,Num_w);
T_max= zeros(1,Num_w);

for i2 = 1:Num_w
    
    v_epsilon=w(i2)*v_t;
    T_max(i2) = fn_find_max(v_epsilon,g);
    
    disp(['w = ',num2str(w(i2)),' max_value = ', num2str(T_max(i2))]);
    
    %     v_k = linspace(-pi,0,num_v);
    %     v_r = zeros(1,num_v);
    %     for i=1:num_v
    %         v_r(i) = fn_m(v_k(i),v_epsilon,g);
    %     end
    %     T_max(i2)=max(v_r);
    
end
hold on;
plot(w,T_max,'*r')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function max_res = fn_find_max(v_epsilon,g)

Num_first_round = 100000;
vk_search = linspace(-pi,0,Num_first_round);
v_res = zeros(1,Num_first_round);
for i=1:Num_first_round
    v_res(i) = fn_m(vk_search(i),v_epsilon,g);
end
[max_res,max_ind] = max(v_res);
center_k = vk_search(max_ind);

Num = 201;
L = pi/2;

while L>1e-8
    
    vk_search = linspace(center_k-L/2,center_k+L/2,Num);
    v_res = zeros(1,Num);
    for i=1:Num
        v_res(i) = fn_m(vk_search(i),v_epsilon,g);
    end
    [max_res,max_ind] = max(v_res);
    center_k = vk_search(max_ind);
    L = L*0.9;
    
end

end


function res = fn_m(k,v_epsilon,g)
Num = length(v_epsilon);
a = 1;

m = eye(2);
for i=1:Num
    m = m*[(2*cos(k*a+1j*g) - v_epsilon(i))*exp(g),-exp(2*g);1,0];
end

res = 2j*sin(k*a+1j*g)*exp(g)/(m(2,1)+m(2,2)*exp(a*(1j*k-2*g/a))...
    -m(1,1)*exp(-1j*k*a)-m(1,2)*exp(-2*g));
res = abs(res)^2;

end


