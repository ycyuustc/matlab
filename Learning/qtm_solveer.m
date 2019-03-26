close(figure(1));
figure(1);

vx = rand(1,4);

vx_delta = rand(1,4)*1e-4;

E = rand();

fn_Fj(vx+vx_delta,E)-fn_Fj(vx,E)
vx_delta*transpose(fn_m(vx,E))

vm = 1:6;
Num = length(vm);
vx = (2*vm-1)/Num/2;

vx = 1j*tan(vx*pi);
beta = 0.1
E = 12/beta;
% E = 99;
vx = fn_descent(vx,E);


v_ba6 = vx/E;

while 1
    vx = fn_descent(vx,E);
    plot([real(vx),1,-1],[imag(vx),0,0],'ro');
    disp(E);
    if E<2
        break;
    end
    E = E-1;
end



disp('finished');

function res = fn_descent(vx,E)

while 1
    
   tm = fn_m(vx,E);
   vFj = fn_Fj(vx,E);
   delta_vx = -tm\transpose(vFj);
   vx = vx + transpose(delta_vx);
   
   tv = fn_Fj(vx,E);
   disp(vx);
   disp(norm(tv));
   
   if norm(tv)<1e-8
       break;
   end
    
end

res = vx;

end


function res = fn_Fj(vx,E)

Num = length(vx);
[mx1,mx2] = meshgrid(vx,vx);
tm = log((mx1-mx2-E)./(mx1-mx2+E));
tm(logical(eye(Num))) = 0;

vx1 = Num*(log((vx-1)./(vx+1))...
    +log((vx+1-E)./(vx-1+E)));
vx2 = sum(tm);
res = vx1 - vx2;

res_real = real(res);
res_imag = imag(res);

res_imag = mod(res_imag+pi,2*pi)-pi;

res = res_real + 1j*res_imag;

end

function res = fn_m(vx,E)
Num = length(vx);
fd = @(x,a) 1./(x-a) - 1./(x+a);
mx2 = ones(Num,1)*vx;
mx1 = transpose(mx2);

tm = fd(mx1-mx2,E);
tm(logical(eye(Num))) = 0;

tv1 = Num*(fd(vx,1)+fd(vx,E-1));
tm1 = diag(tv1);
tv2 = transpose(sum(tm,2));
tm2 = diag(tv2);

tm3 = tm;

res = tm1 - tm2 + tm3;

end