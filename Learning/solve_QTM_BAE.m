% function solve_QTM_BAE

Num = 4;

v_initial = -(Num-1):2:(Num-1);

vx = v_initial*0.5002;
tau = 0.00j;
fn_f1 = @(x,n) power((x+tau)./(x-tau).*(x-1-tau)./(x+1+tau),n);

disp(fn_f2(vx,2));

disp(vx);
while 1
    
    vx_pre = vx;
    
    for k=1:Num
        
        if k==1
            left_bound = -100;
            right_bound = vx(k+1)-1;
        else
            if k==Num
                left_bound = vx(k-1)+1;
                right_bound = 100;
            else
                left_bound = vx(k-1)+1;
                right_bound = vx(k+1)-1;
            end
        end
        
        mystep = 0.0001;
        while 1
            
            tvx1 = vx;
            tvx2 = vx;
            
            while 1
                test_left = vx(k) - mystep;
                test_right = vx(k) + mystep;
                
                if test_left <= left_bound
                    mystep = mystep*0.5;
                end
                
                if test_right >= right_bound
                    mystep = mystep*0.5;
                end
                
                if test_left>left_bound && test_right<right_bound
                    break;
                end
            
            end
            
            tvx1(k) = test_left;
            tvx2(k) = test_right;
            
            dis_0 = fn_judge(vx, Num);
            dis_1 = fn_judge(tvx1,Num);
            dis_2 = fn_judge(tvx2,Num);
            
            if dis_1<dis_0
                vx = tvx1;
%                 disp('left_move');
            end
            
            if dis_2<dis_0
                vx = tvx2;
%                 disp('right_move');
            end
            
%             disp(vx);
            
            if dis_1>=dis_0 && dis_2>=dis_0
                mystep = mystep*0.5;
%                 disp(dis_0);
            end
            
            
            if mystep<1e-8
                break;
            end
            
        end
        
    end
    
    disp(vx);
    
    if norm(vx_pre - vx)<1e-8
        
        break;
        
    end
    
    
end

for i=1:Num
    disp(fn_Fj(vx,i,Num));
end

disp(fn_f1(vx(3),Num));
disp(fn_f2(vx,3));


% end

function res = fn_f2(vx, ind)
x = vx(ind);
tv1 = x-1-vx;
tv2 = x+1-vx;
tv = tv1./tv2;
res = -prod(tv);

end

function res = fn_Fj(vx,ind,Num)

tau = 0.00j;
fn_f1 = @(x,n) power((x+tau)./(x-tau).*(x-1-tau)./(x+1+tau),n);
lhs = fn_f1(vx(ind),Num);
rhs = fn_f2(vx,ind);
res = abs(lhs - rhs);

end

function [res,v_lhs,v_rhs] = fn_judge(vx,Num)
tau = 0;
fn_f1 = @(x,n) power((x+tau)./(x-tau).*(x-1-tau)./(x+1+tau),n);
[mx,my] = meshgrid(vx,vx);
tm = (mx-my-1)./(mx-my+1);
tm(logical(eye(Num))) = 1;
v_lhs = fn_f1(vx,Num);
v_rhs = prod(tm);
res = norm(v_lhs./v_rhs-1);
% disp(v_lhs./v_rhs);

end



