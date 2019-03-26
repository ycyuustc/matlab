Num = 20;

Num_solve = 6;
Num_fixed = Num - Num_solve;

vx = -(Num_fixed-1):2:(Num_fixed-1);
vx = vx/2;

v_right = 2*(1:(Num_solve/2))+Num_fixed/2-1/2;
v_left = -v_right(end:-1:1);
vx = [v_left,vx,v_right];

v_ind_solve = find(abs(vx)>Num_fixed/2);

while 1
    
    tm = fn_m(vx);
    vFj = fn_Fj(vx);
    
    tm = tm(v_ind_solve,v_ind_solve);
    vFj = vFj(v_ind_solve);
    
    delta_vx = -tm\transpose(vFj);
    vx(v_ind_solve) = vx(v_ind_solve) + transpose(delta_vx);
    %         disp(vx);
    tv = fn_Fj(vx);
    tv = tv(v_ind_solve);
    disp(norm(tv));
    disp(vx);
    
    if norm(tv)<1e-8
        break;
    end
        
end



function res = fn_Fj(vx)
Num = length(vx);
[mx1,mx2] = meshgrid(vx,vx);

f2 = @(x) log((x-1)./(x+1));

tm = f2(mx1-mx2);
tm(logical(eye(Num))) = 0;
vx1 = Num*f2(vx);
vx2 = sum(tm);
res = vx1-vx2;

end

function res = fn_m(vx)

Num = length(vx);
d2 = @(x) 1./(x-1) - 1./(x+1);

mx2 = ones(Num,1)*vx;
mx1 = transpose(mx2);

tm = d2(mx1-mx2);
tm(logical(eye(Num))) = 0;

tv1 = Num*d2(vx);
tm1 = diag(tv1);
tv2 = transpose(sum(tm,2));
tm2 = -diag(tv2);

tm3 = tm;

res = tm1 + tm2 + tm3;

end
