% function solve_QTM_BAE
tau = 0;
fn_f1 = @(x,n) power((x+tau)./(x-tau).*(x-1-tau)./(x+1+tau),n);
Num = 4;

v_initial = -(Num-1):2:(Num-1);

vx = v_initial*0.6;
fn_judge(vx,Num)

ind0 = 3;
[ind1,ind2] = fn_closest(vx,ind0)

tvx = vx;
tvx([ind1,ind2,ind0]) = [];

a = fn_rhs(vx(ind0),tvx)



ind0 = 1;
while 1
    
    [ind1,ind2] = fn_closest(vx,ind0); 
    tvx = vx;
    tvx([ind1,ind2,ind0]) = [];
    
    x0 = vx(ind0);
    x1 = vx(ind1);
    x2 = vx(ind2);
    a = fn_f1(x0,Num)/fn_rhs(x0,tvx);
    
    [try_1,try_2] = fn_solve_step(x1,x2,a);
    
    if abs(try_1 - x0)<abs(try_2 - x0)
        vx(ind0) = try_1;
    else
        vx(ind0) = try_2;
    end
    
    disp(vx);
    
    delta = fn_judge(vx,Num);
    disp(delta);
    if delta<1e-8
        break;
    end
    
    ind0 = ind0 + 1;
    if ind0>Num
        ind0 = ind0 - Num;
    end
    
end




% end

function [ind1,ind2] = fn_closest(vx,ind)

x = vx(ind);
v_delta = abs(vx-x);
[~,v_sort] = sort(v_delta);
ind1 = v_sort(2);
ind2 = v_sort(3);

end

function res = fn_rhs(x,vy)

res = (x-vy-1)./(x-vy+1);
res = prod(res);

end

function res = fn_judge(vx,Num)
tau = 0;
fn_f1 = @(x,n) power((x+tau)./(x-tau).*(x-1-tau)./(x+1+tau),n);
[mx,my] = meshgrid(vx,vx);
tm = (mx-my-1)./(mx-my+1);
tm(logical(eye(Num))) = 1;
v_lhs = fn_f1(vx,Num);
v_rhs = prod(tm);
res = norm(v_lhs./v_rhs-1);
disp(v_lhs./v_rhs);

end



function res = fn_f2(vx, ind)
x = vx(ind);
tv1 = x-1-vx;
tv2 = x+1-vx;
tv = tv1./tv2;
res = -prod(tv);

end

function res = fn_Fj(vx,ind,Num)

tau = 0.001j;
fn_f1 = @(x,n) power((x+tau)./(x-tau).*(x-1-tau)./(x+1+tau),n);
lhs = fn_f1(vx(ind),Num);
rhs = fn_f2(vx,ind);
res = abs(lhs - rhs);

end

function [res1,res2] = fn_solve_step(x1,x2,a)

Delta = 16*a+x1^2-2*a*x1^2+a^2*x1^2-2*x1*x2+4*a*x1*x2-2*a^2*x1*x2...
    +x2^2-2*a*x2^2+a^2*x2^2;

Dsq = sqrt(Delta);
res1 = (2+2*a+x1-a*x1+x2-a*x2+Dsq)/(2*(1-a));
res2 = (-2-2*a-x1+a*x1-x2+a*x2+Dsq)/(2*(-1+a));

end



