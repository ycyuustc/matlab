Num = 16;
vm = (0:(Num/2))+1e-4;
delta_b = 0.5*j;
b = (1/2+delta_b);

v_dd = -exp(4*1j*pi*vm/Num)-b^2*(-1+exp(4*1j*pi*vm/Num)).^2;
v_dd = sqrt(v_dd);
vx1 = -(1j+1j*exp(4*1j*pi*vm/Num)+2*v_dd)./(2*(-1+exp(4*1j*vm*pi/Num)));
vx2 = -(1j+1j*exp(4*1j*pi*vm/Num)-2*v_dd)./(2*(-1+exp(4*1j*vm*pi/Num)));
disp(vx1);
disp(vx2);

vx = [vx1,vx2];

close(figure(1));
close(figure(2));
figure(1);
hold on;
plot(real(vx1),imag(vx1),'bo');
plot(real(vx2),imag(vx2),'bo');
plot(0,1,'ro');
plot(0,-1,'ro');
epsilon = 1e6;
gamma = 1e-6;
fn_Fj(vx,b,epsilon)

tm = fn_m(vx,b,epsilon);

v_delta_x = rand(1,Num)*1e-5;

fn_Fj(vx+v_delta_x,b,epsilon) - fn_Fj(vx,b,epsilon)
v_delta_x*transpose(fn_m(vx,b,epsilon))

vx0 = vx;

while 1
    while 1
        
        tm = fn_m(vx,b,epsilon);
        vFj = fn_Fj(vx,b,epsilon);
        
        delta_vx = -tm\transpose(vFj);
        vx = vx + transpose(delta_vx);
%         disp(vx);
        disp(norm(fn_Fj(vx,b,epsilon)));
        
        if norm(fn_Fj(vx,b,epsilon))<1e-8
            break;
        end
        
    end
    
%     figure(1);
%     hold on;
%     plot(real(vx),imag(vx),'g*');
    
    figure(2);
    plot(real(vx),imag(vx),'ro');
    
    pause(0.001);
    
    
    gamma = gamma + 1e-3;
    epsilon = 1/gamma;
    disp(epsilon);
    
    if epsilon<1.18
        disp('look out!');
    end
    
end

function res = fn_Fj(vx,b,epsilon)
Num = length(vx);
[mx1,mx2] = meshgrid(vx,vx);

f1 = @(x) log(((x-1j/2).^2+b^2)./((x+1j/2).^2+b^2));
f2 = @(x) log((x-1j*epsilon)./(x+1j*epsilon));

tm = f2(mx1-mx2);
tm(logical(eye(Num))) = 0;
vx1 = Num/2*f1(vx);
vx2 = sum(tm);
res = vx1-vx2;

res_real = real(res);
res_imag = mod(imag(res)+pi,2*pi) - pi;
res = res_real + 1j*res_imag;
end

function res = fn_m(vx,b,epsilon)

Num = length(vx);
d2 = @(x) 1./(x-1j*epsilon) - 1./(x+1j*epsilon);
d1 = @(x) -8j*(-1+4*b^2-4*x.^2)./(16*b^4+8*b^2*(-1+4*x.^2)+(1+4*x.^2).^2);

mx2 = ones(Num,1)*vx;
mx1 = transpose(mx2);

tm = d2(mx1-mx2);
tm(logical(eye(Num))) = 0;

tv1 = Num/2*d1(vx);
tm1 = diag(tv1);
tv2 = transpose(sum(tm,2));
tm2 = -diag(tv2);

tm3 = tm;

res = tm1 + tm2 + tm3;

end
