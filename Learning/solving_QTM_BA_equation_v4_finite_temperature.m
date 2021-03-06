f = @(x) log(abs((x-1)./(x+1)));
Num = 20;
epsilon = 0;
Num_solve = 6;
Num_fixed = Num - Num_solve;

m_d2x = triu(ones(Num));
m_x2d = eye(Num) - diag(ones(1,Num-1),1);

m_d2x_b = diag([-ones(1,Num-1),1])*(triu(ones(Num))');
m_x2d_b = diag(ones(1,Num-1),-1)+diag([-ones(1,Num-1),1]);
mM = 2*eye(Num-1) - diag(ones(1,Num-2),1) - diag(ones(1,Num-2),-1);

vx = -(Num_fixed-1):2:(Num_fixed-1);
vx = vx/2;

v_right = 2*(1:(Num_solve/2))+Num_fixed/2-1/2;
v_left = -v_right(end:-1:1);
vx = [v_left,vx,v_right];

v_ind_solve = find(abs(vx)>Num_fixed/2);

vx = fn_descent(vx,epsilon,v_ind_solve);
vx0 = vx;
vL = zeros(1,Num);

close(figure(1));
figure(1);
hold on;

while 1

while 1
    vx_pre = vx;
    vL_pre = vL;
    for p=1:Num
        
        tv = vx;
        if p>=2 && p<=Num-1
            tv((p-1):(p+1))=[];
        else
            if p==1
                tv(1:2) = [];
            end
            if p==Num
                tv((Num-1):Num) = [];
            end
        end
        
        x = vx(p);
        vL(p) = Num*f(x)-sum(f(x-tv));
        
    end
    
    vL = 1/2*(vL - vL(end:-1:1));
    disp(' vL = ');
    disp(vL);
    disp(vL'-vL_pre');
    
    vN = vL(2:end) - vL(1:end-1);
    vN = transpose(vN);
    
    vD = -mM\vN;
    
    disp('vD = ');
    disp(vD);
    
    vd = 1./tanh(vD/2);
    
    vd = 1/2*(vd + vd(end:-1:1));
       
    disp('vd = ');
    disp(vd);
    
    tv_ind = find(vd<1);
    vd(tv_ind) = 2;
    
    x0 = 0 - sum(vd)/2;
    
    vx = [x0,vd']*m_d2x;
    disp(vx);
    
    vx = fn_descent(vx,epsilon,fn_tobesolve(vx));
    
    disp('vx');
    disp(vx);
    
    disp(epsilon);
    disp('norm');
    disp(norm(vx_pre-vx));
    
    if norm(vx_pre-vx)<1e-8
        break;
    end
    
end

disp(fn_Fj(vx,epsilon));

plot(vx,ones(1,Num)*epsilon,'ro');

epsilon = epsilon + 1e-3;

if epsilon>0.1
    break;
end

end


% =================== END of MAIN SCRIPTS =====================

function res = fn_tobesolve(vx)
Num = length(vx);
mx2 = ones(Num,1)*vx;
mx1 = transpose(mx2);

tm = abs(abs(mx1-mx2)-1);

tv = min(tm);

non_singularity_index = find(tv>1e-2);
res = non_singularity_index;

end

function res = fn_descent(vx,epsilon,v_ind_solve)

vFj = fn_Fj(vx,epsilon);
tv = vFj(v_ind_solve);
norm_pre = norm(tv);

while 1
    
    step = 1;
    while 1
    vx = (vx - vx(end:-1:1))/2;
    vx0 = vx;
    
    tm = fn_m(vx,epsilon);
    vFj = fn_Fj(vx,epsilon);
    
    tm = tm(v_ind_solve,v_ind_solve);
    vFj = vFj(v_ind_solve);
    
    vFj = (vFj - vFj(end:-1:1))/2;
    
    delta_vx = -tm\transpose(vFj);
    vx(v_ind_solve) = vx(v_ind_solve) + step*transpose(delta_vx);
    
    tv = fn_Fj(vx,epsilon);
    tv = tv(v_ind_solve);
    norm_now = norm(tv);
    
    disp(norm_now);
    disp(vx);
    
    if(norm_now>norm_pre)
        if step>1e-6
            step = step/2;
            vx = vx0;
        else
            disp('hehe');
            error('error');
        end
    else
        norm_pre = norm_now;
        break;
    end
    
    end
    
    if norm(tv)<1e-8
        break;
    end
end

res = vx;

end


function res = fn_Fj(vx,epsilon)

Num = length(vx);
[mx1,mx2] = meshgrid(vx,vx);
tm = log(abs((mx1-mx2-1)./(mx1-mx2+1)));
tm(logical(eye(Num))) = 0;

vx1 = Num*(log(abs((vx-1+epsilon)./(vx+1-epsilon)))...
    -log(abs((vx-epsilon)./(vx+epsilon))));
vx2 = sum(tm);
res = vx1 - vx2;

end

function res = fn_m(vx,epsilon)
Num = length(vx);
fd = @(x,a) 1./(x-a) - 1./(x+a);
mx2 = ones(Num,1)*vx;
mx1 = transpose(mx2);

tm = fd(mx1-mx2,1);
tm(logical(eye(Num))) = 0;

tv1 = Num*(fd(vx,1-epsilon)+fd(vx,epsilon));
tm1 = diag(tv1);
tv2 = transpose(sum(tm,2));
tm2 = -diag(tv2);

tm3 = tm;

res = tm1 + tm2 + tm3;

end

