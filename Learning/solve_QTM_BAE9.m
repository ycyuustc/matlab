% 
% 
% close(figure(1));
% figure(1);
% hold on;
% close(figure(2));
% figure(2);
% hold on;
% close(figure(3));
% figure(3);
% hold on;

tau = 0;
fn_f1 = @(x,n) power((x+tau)./(x-tau).*(x-1-tau)./(x+1+tau),n);
Num = 6;

v_initial = -(Num-1):2:(Num-1);

vx = v_initial*0.6;

x1 = vx(1);
x2 = vx(2);
x3 = vx(3);
x4 = vx(4);

tv = fn_force(vx,Num);

delta_T = 4e-3;
mu1 = 10.0;
mu2 = 0.05;
mu3 = 0.01;
mass = 4.0;
distance = 1.0;

vv = zeros(1,Num);
disp(vx);
vx_pre = vx;
y = 0;
while 1
    
    [jm,vx] = fn_judge_matrix(vx,vx_pre,distance);
    tv = sum(jm);
    tv = 1./tv;
    re_jm = diag(tv)*jm;
    mark_jm = 1 - jm;
       
    [lhs,fm] = fn_force_matrix(vx,Num);
    lhs = lhs*re_jm;
    fm = fm.*mark_jm;
    rhs = sum(fm);
    rhs = rhs*re_jm;
    vv = vv*re_jm;
    
    va = (lhs + rhs)/mass;
    
%     va = fn_force(vx,Num)/mass;
    va = va ... %- mu1*abs(sqrt(vv)).*sign(vv) ...
         - 0.2*(exp(abs(vv)/mu2)-1).*sign(vv);
    %     -mu3.*sign(vv);
    
    vx_pre = vx;
    vx = vx + delta_T*vv;
    vv = vv + delta_T*va;
    
    disp(vx);
    
%     disp(['x:           ',num2str(vx)]);
%     disp(['speed:       ',num2str(vv)]);
%     disp(['accelerate:  ',num2str(va)]);
    
%     figure(1);
%     plot(vx,y*ones(1,Num),'b*');
%     figure(2);
%     plot(vv,y*ones(1,Num),'g*');
%     figure(3);
%     plot(va,y*ones(1,Num),'r*');
%     figure(1);
%     plot(vx(1),y,'b*');
%     plot(vx(2),y,'g*');
%     plot(vx(3),y,'r*');
%     plot(vx(4),y,'y*');
%     plot(vx(5),y,'b*');
%     plot(vx(6),y,'g*');
%     figure(2);
%     plot(vv(1),y,'b*');
%     plot(vv(2),y,'g*');
%     plot(vv(3),y,'r*');
%     plot(vv(4),y,'y*');
%     plot(vv(5),y,'b*');
%     plot(vv(6),y,'g*');
%     figure(3);
%     plot(va(1),y,'b*');
%     plot(va(2),y,'g*');
%     plot(va(3),y,'r*');
%     plot(va(4),y,'y*');
%     plot(va(5),y,'b*');
%     plot(va(6),y,'g*');
    
    
    y = y + delta_T;
    
%     pause(0.01);
        
end



% end


function [matrix,tvx2] = fn_judge_matrix(tvx,tvx_pre,distance)

Num = length(tvx);
tvx2 = fn_normalize(tvx,tvx_pre,distance);

[mx,my] = meshgrid(tvx2,tvx2);
tm = zeros(Num,Num);
matrix = zeros(Num,Num);
tm(abs(mx-my+distance)<1e-6) = 1;

for k=1:Num
    
    matrix = matrix+tm^k;
    
end

matrix = matrix + matrix';

matrix(logical(eye(Num))) = 1;

end

function vx = fn_normalize(vx,vx_pre,distance)

Num = length(vx);
[mx,my] = meshgrid(vx,vx);
[mx_pre,my_pre] = meshgrid(vx_pre,vx_pre);

m_delta = -mx+my;
m_delta_pre = -mx_pre + my_pre;

m1 = sign(m_delta + distance);
m2 = sign(m_delta_pre + distance);
m3 = m1.*m2;

m4 = zeros(Num,Num);
m5 = ones(Num,Num);

m4(m3<=0) = 1;
m5(m2==0) = 0;

m6 = m4.*m5;

v6 = sum(m6,2);

for k=1:Num
   
    if v6(k)>0
       
        tv = vx(logical(m6(k,:)));
        pnt = min(tv);
        DistanceNow = pnt - vx(k);
        NeedToMove = distance - DistanceNow;
        
        vx(vx<vx(k)) = vx(vx<vx(k)) - NeedToMove;
        vx(k) = vx(k) - NeedToMove;
        
    end
    
end

end


function [matrix,tvx2] = fn_judge_matrix_old(tvx,distance)
vx = tvx;
Num = length(tvx);
tvx2 = tvx;
for i=1:(Num-1)
    
    delta = tvx(i+1) - tvx(i);
    if delta<=distance + 1e-6
        tvx(i+1:end) = tvx(i+1:end) - delta;
        tvx2(i+1:end) = tvx2(i+1:end) + (distance - delta)/2.0;
        tvx2(1:i) = tvx2(1:i) - (distance - delta)/2.0;
    end
    
end

[mx,my] = meshgrid(vx,vx);
matrix = zeros(Num,Num);
delta_m = abs(mx-my);
matrix(abs(floor(delta_m+0.5)-delta_m)<1e-3) = 1;

end

function [lhs,mr] = fn_force_matrix(vx,Num)

fn_f1=@(x) log(abs((x-1)./(x+1)));
fn_f3=@(x) log(abs((x-1)./(x+1)));

[mx,my] = meshgrid(vx,vx);

lhs = Num*fn_f1(vx);
mr = -fn_f3(mx-my);
mr(logical(eye(Num))) = 0;

mr(mr==Inf) = 0;
mr(mr==-Inf) = 0;
% rhs = sum(mr);

% res = lhs + rhs;

end



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



function [res,lhs,mr] = fn_force(vx,Num)

fn_f1=@(x) log(abs((x-1)./(x+1)));
fn_f3=@(x) log(abs((x-1)./(x+1))).*(sign(x.^2-1)+1)/2 ...
     +(sign(1-x.^2)+1)/2.*sign(-x)*20;

[mx,my] = meshgrid(vx,vx);

lhs = Num*fn_f1(vx);
mr = -fn_f3(mx-my);
mr(logical(eye(Num))) = 0;
rhs = sum(mr);

res = lhs + rhs;

end



