
beta = 1e-1;
M = 12;
N = M/2;
h = 1e-1;

v_mu = fn_solve(beta,M,h);

epsilon = beta/M;

tv = prod((v_mu*epsilon+1)./v_mu)*exp(beta*h/2)...
    +prod((1-epsilon*v_mu)./v_mu)*exp(-beta*h/2);

tv = tv/(1-epsilon)^N;
disp(tv);



disp('finished');

function res = fn_descent(vx,E,beta_h)

while 1
    
   tm = fn_m(vx,E);
   vFj = fn_Fj(vx,E,beta_h);
   delta_vx = -tm\transpose(vFj);
   vx = vx + transpose(delta_vx);
   
   tv = fn_Fj(vx,E,beta_h);
   disp(vx);
   disp(norm(tv));
   
   if norm(tv)<1e-10
       break;
   end
    
end

res = vx;

end


function res = fn_solve(beta,M,h)

Num = M/2;
E = M/beta;

vm = 1:Num;

vx = (2*vm-1)/Num/2;

vx = 1j*tan(vx*pi);

res = fn_descent(vx,E,beta*h);

end

function res = fn_Fj(vx,E,beta_h)

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

res = res - beta_h;

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