% Num = 8;
% tv1 = -20.5:20.5;
% Num = length(tv1);
% vx0 = linspace(-0.1,0.1,Num)+tv1;
% vd0 = rand(1,Num)+1;
f = @(x) log(abs((x-1)./(x+1)));

m_d2x = triu(ones(Num));
m_x2d = eye(Num) - diag(ones(1,Num-1),1);

m_d2x_b = diag([-ones(1,Num-1),1])*(triu(ones(Num))');
m_x2d_b = diag(ones(1,Num-1),-1)+diag([-ones(1,Num-1),1]);
mM = 2*eye(Num-1) - diag(ones(1,Num-2),1) - diag(ones(1,Num-2),-1);

% vx0 = [-1.502,-0.501,0.5,1.504];


% while 1
%     
%    tm = fn_m(vx);
%    vFj = fn_Fj(vx);
%    
%    delta_vx = -tm\transpose(vFj);
%    vx = vx+transpose(delta_vx);
%    disp(vx);
%    disp(norm(fn_Fj(vx)));
%     
% end

Num = 10;

Num_solve = 8;
Num_fixed = Num - Num_solve;

vx = -(Num_fixed-1):2:(Num_fixed-1);
vx = vx/2;

v_right = 2*(1:(Num_solve/2))+Num_fixed/2-1/2;
v_left = -v_right(end:-1:1);
vx = [v_left,vx,v_right];

v_ind_solve = find(abs(vx)>Num_fixed/2);

vx = fn_descent(vx,v_ind_solve);
vx0 = vx;
vL = zeros(1,Num);
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
    
    vx = fn_descent(vx,fn_tobesolve(vx));
    
    disp('vx');
    disp(vx);
    
    disp('norm');
    disp(norm(vx_pre-vx));
    
    if norm(vx_pre-vx)<1e-8
        break;
    end
    
end

disp(fn_Fj(vx));

% while 1
%     
%     vx_total_pre = vx;
%     vd = vx*m_x2d;
%     for point = 1:(Num-1)
%         while 1
%             vx_pre = vx;
%             x = vx(point);
%             tv = vx;
%             tv(point) = [];
%             tv(point) = [];
%             LHS = 4*Num/2*f(x);
%             RHS = sum(f(x-tv));
%             d = 1 + 2/(exp(LHS-RHS)-1);
%             vd(point+1) = d; % d_point update !
%             vx = vd*m_d2x;
%             disp(point);
%             disp(vx);
%             disp(vd);
%             disp(fn_Fj(vx));
%             
%             if norm(vx_pre-vx)<1e-8
%                 break;
%             end
%         end
%         
%     end
%     
%     
%     vd_b = vx*m_x2d_b;
%     
%     for point = Num:-1:2
%         while 1
%             vx_pre = vx;
%             x = vx(point);
%             tv = vx;
%             tv(point) = [];
%             tv(point-1) = [];
%             LHS = 4*Num/2*f(x);
%             RHS = sum(f(x-tv));
%             d = 1 + 2/(exp(-LHS+RHS)-1);
%             vd_b(point-1) = d; % d_(point-1) update
%             vx = vd_b*m_d2x_b;
%             disp(point);
%             disp(vx);
%             disp(vd_b);
%             disp(fn_Fj(vx));
%             
%             if norm(vx_pre-vx)<1e-8
%                 break;
%             end
%         end
%         
%     end
%     
%     disp(vx);
%     disp(norm(vx_total_pre-vx));
%     
%     if norm(vx_total_pre-vx)<1e-8
%         break;
%     end
%     
% end
% 
% 
% step = 0.1;
% while step>1e-8
%     vx_left = vx - step*ones(1,Num);
%     vx_right = vx + step*ones(1,Num);
%     
%     value_left = sum(abs(fn_Fj(vx_left)));
%     value_right = sum(abs(fn_Fj(vx_right)));
%     value_center = sum(abs(fn_Fj(vx)));
%     
%     if value_center<min(value_left,value_right)
%         step = step/2;
%         disp(step);
%     else
%         if value_left<value_right
%             vx = vx_left;
%         else
%             vx = vx_right;
%         end
%     end
%     
%     disp(vx);
%     disp(sum(abs(fn_Fj(vx_right))));
%     disp(vx(2:end)-vx(1:end-1));
%     
% end

function res = fn_tobesolve(vx)
Num = length(vx);
mx2 = ones(Num,1)*vx;
mx1 = transpose(mx2);

tm = abs(abs(mx1-mx2)-1);

tv = min(tm);

non_singularity_index = find(tv>1e-3);
res = non_singularity_index;

end

function res = fn_descent(vx,v_ind_solve)

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

res = vx;

end


function res = fn_Fj(vx)

Num = length(vx);
[mx1,mx2] = meshgrid(vx,vx);
tm = log(abs((mx1-mx2-1)./(mx1-mx2+1)));
tm(logical(eye(Num))) = 0;

vx1 = Num*log(abs((vx-1)./(vx+1)));
vx2 = sum(tm);
res = vx1 - vx2;

end

function res = fn_m(vx)
Num = length(vx);
fd = @(x) 1./(x-1) - 1./(x+1);
mx2 = ones(Num,1)*vx;
mx1 = transpose(mx2);

tm = fd(mx1-mx2);
tm(logical(eye(Num))) = 0;

tv1 = Num/2*fd(vx);
tm1 = diag(tv1);
tv2 = transpose(sum(tm,2));
tm2 = -diag(tv2);

tm3 = tm;

res = tm1 + tm2 + tm3;

end

function res = fn_Fj_dx0(vx)

Num = length(vx);
res = Num/2*(1./(abs(vx-1)) - 1./(abs(vx+1)));

end