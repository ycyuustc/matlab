vx = rand(1,3);
Lambda = rand();
c0 = rand();
c1 = rand();
vm = [0,-1,1]; 

vx = rand(1,3);
vx_delta = rand(1,3)*1e-4;

vx = 2*pi*vm;

% fn_Fj(vx+vx_delta,Lambda,c0,c1,vm)-fn_Fj(vx,Lambda,c0,c1,vm)
% vx_delta*transpose(fn_m(vx,Lambda,c0,c1))
c1 = -10;
Num_step = 100;
c0 = 0.01;
v_c1 = linspace(0,c1,Num_step);

for i=2:Num_step
   c1_now = v_c1(i);
   vx = fn_steepest(vx,0,c0,c1_now,vm); 
end

fn_Fj(vx,0,c0,c1,vm)

Num_step = 100;
v_c0 = linspace(0,c0,Num_step);
v_c1 = linspace(0,c1,Num_step);

for i=2:Num_step
   c0_now = v_c0(i);
   vx = fn_steepest(vx,K,c0_now,0,vm);
   vk = vx(1:3);
   plot(real(vk),imag(vk),'ro');
   
end


disp(fn_Fj(vx,K,c0,0,vm));
for i=2:Num_step
   c1_now = v_c1(i);
   vx = fn_steepest2(vx,K,c0,c1_now,vm);
   vk = vx(1:3);
   plot(real(vk),imag(vk),'ro');
   
end



function res = fn_steepest(vx0,Lambda,c0,c1,vm)

vx = vx0;
while 1
   
    tm = fn_m(vx,Lambda,c0,c1);
    vFj = fn_Fj(vx,Lambda,c0,c1,vm);
    delta_vx = -tm\transpose(vFj);
    vx = vx + transpose(delta_vx);
    
    norm_now = norm(fn_Fj(vx,Lambda,c0,c1,vm));
    disp(vx);
    disp(norm_now);
    
    if norm_now<1e-8
        break;
    end
    
    
end

res = vx;

end

function res = fn_steepest2(vx0,K,c0,c1,vm)

vx = vx0;
while 1
   
    tm = fn_m2(vx,K,c0,c1);
    vFj = fn_Fj2(vx,K,c0,c1,vm);
    delta_vx = -tm\transpose(vFj);
    vx = vx + transpose(delta_vx);
    
    norm_now = norm(fn_Fj2(vx,K,c0,c1,vm));
    disp(vx);
    disp(norm_now);
    
    if norm_now<1e-8
        break;
    end
    
    
end

res = vx;

end





function res = fn_Fj(vx,Lambda,c0,c1,vm)

c = @(x2,x3) c0 + c1./(1+exp(-x2-x3));
f = @(x,y) log((x-Lambda+1j/2*y)/(x-Lambda-1j/2*y));
vFj = zeros(1,3);
x1 = vx(1); x2 = vx(2); x3 = vx(3);

vFj(1) = 1j*x1 - f(x1,c(x2,x3)) - vm(1)*1j*2*pi;
vFj(2) = 1j*x2 - f(x2,c(x3,x1)) - vm(2)*1j*2*pi;
vFj(3) = 1j*x3 - f(x3,c(x1,x2)) - vm(3)*1j*2*pi;

res = vFj;

end

function res = fn_m(vx,Lambda,c0,c1)
c = @(x2,x3) c0 + c1./(1+exp(-x2-x3));
dc = @(x2,x3) c1*exp(-x2-x3)/(1+exp(-x2-x3))^2;

m = zeros(3,3);
x1 = vx(1); x2 = vx(2); x3 = vx(3);

m(1,1) = 1j - 1/(x1-Lambda+1j/2*c(x2,x3)) + 1/(x1-Lambda-1j/2*c(x2,x3));
m(1,2) = -1j/2*dc(x2,x3)/(x1-Lambda+1j/2*c(x2,x3)) - 1j/2*dc(x2,x3)/(x1-Lambda-1j/2*c(x2,x3));
m(1,3) = -1j/2*dc(x2,x3)/(x1-Lambda+1j/2*c(x2,x3)) - 1j/2*dc(x2,x3)/(x1-Lambda-1j/2*c(x2,x3));

m(2,2) = 1j - 1/(x2-Lambda+1j/2*c(x3,x1)) + 1/(x2-Lambda-1j/2*c(x3,x1));
m(2,1) = -1j/2*dc(x3,x1)/(x2-Lambda+1j/2*c(x3,x1)) - 1j/2*dc(x3,x1)/(x2-Lambda-1j/2*c(x3,x1));
m(2,3) = -1j/2*dc(x3,x1)/(x2-Lambda+1j/2*c(x3,x1)) - 1j/2*dc(x3,x1)/(x2-Lambda-1j/2*c(x3,x1));

m(3,3) = 1j - 1/(x3-Lambda+1j/2*c(x1,x2)) + 1/(x3-Lambda-1j/2*c(x1,x2));
m(3,1) = -1j/2*dc(x1,x2)/(x3-Lambda+1j/2*c(x1,x2)) - 1j/2*dc(x1,x2)/(x3-Lambda-1j/2*c(x1,x2));
m(3,2) = -1j/2*dc(x1,x2)/(x3-Lambda+1j/2*c(x1,x2)) - 1j/2*dc(x1,x2)/(x3-Lambda-1j/2*c(x1,x2));

res = m;

end



function res = fn_Fj2(vx,K,c0,c1,vm)

vx(4) = 1/vx(4);

c = @(x) c0 + c1./(1+exp(x-K));
f = @(x1,x2) log((x1+1j/2*c(x2))./(x1-1j/2*c(x2)));

vFj = zeros(1,4);

vFj(1) = 1j*vx(1) - f(vx(1)-vx(4),vx(3)) - 1j*vm(1)*2*pi;
vFj(2) = 1j*vx(2) - f(vx(2)-vx(4),vx(1)) - 1j*vm(2)*2*pi;
vFj(3) = 1j*vx(3) - f(vx(3)-vx(4),vx(2)) - 1j*vm(3)*2*pi;
vFj(4) = vx(1)+vx(2)+vx(3)-K;

res = vFj;

end

function res = fn_m2(vx,K,c0,c1)
c = @(x) c0 + c1./(1+exp(x-K));
dc = @(x) -c1.*exp(x-K)./((1+exp(x-K)).^2);
f = @(x,y) 1./(x+y);
m = zeros(4,4);

vx(4) = 1/vx(4);

x1 = vx(1); x2 = vx(2); x3 = vx(3); Lambda = vx(4);

m(1,1) = 1j - f(x1-Lambda,1j/2*c(x3)) + f(x1-Lambda,-1j/2*c(x3));
m(1,3) = (-f(x1-Lambda,1j/2*c(x3)) - f(x1-Lambda,-1j/2*c(x3)))*1j/2*dc(x3);
m(1,4) = (f(x1-Lambda,1j/2*c(x3)) - f(x1-Lambda,-1j/2*c(x3)))*(-Lambda^2);

m(2,2) = 1j - f(x2-Lambda,1j/2*c(x1)) + f(x2-Lambda,-1j/2*c(x1));
m(2,1) = (-f(x2-Lambda,1j/2*c(x1)) - f(x2-Lambda,-1j/2*c(x1)))*1j/2*dc(x1);
m(2,4) = (f(x2-Lambda,1j/2*c(x1)) - f(x2-Lambda,-1j/2*c(x1)))*(-Lambda^2);

m(3,3) = 1j - f(x3-Lambda,1j/2*c(x2)) + f(x3-Lambda,-1j/2*c(x2));
m(3,2) = (-f(x3-Lambda,1j/2*c(x2)) - f(x3-Lambda,-1j/2*c(x2)))*1j/2*dc(x2);
m(3,4) = (f(x3-Lambda,1j/2*c(x2)) - f(x3-Lambda,-1j/2*c(x2)))*(-Lambda^2);

m(4,:) = [1,1,1,0];

res = m;

end