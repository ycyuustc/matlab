v_epsilon = rand(1,40);
v_k = linspace(-pi,0,1000);
v_r = zeros(1,1000);
for i=1:1000
   v_r(i) = fn_m(v_k(i),v_epsilon); 
end
close(figure(10));figure(10);
plot(v_k,v_r,'r*');

function res = fn_m(k,v_epsilon)
Num = length(v_epsilon);
a = 1;
g = -0.3;
m = eye(2);
for i=1:Num
   m = m*[(2*cos(k*a)+1j*g - v_epsilon(i))*exp(g),-exp(2*g);1,0]; 
end

res = 2j*sin(k*a+1j*g)*exp(g)/(m(2,1)+m(2,2)*exp(a*(1j*k-2*g/a))...
    -m(1,1)*exp(-1j*k*a)-m(1,2)*exp(-2*g));
res = abs(res)^2;

end