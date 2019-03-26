Num = 4;
% vd0 = rand(1,Num)+1;
f = @(x) log(abs((x-1)./(x+1)));

m_d2x = triu(ones(Num));
m_x2d = eye(Num) - diag(ones(1,Num-1),1);

m_d2x_b = diag([-ones(1,Num-1),1])*(triu(ones(Num))');
m_x2d_b = diag(ones(1,Num-1),-1)+diag([-ones(1,Num-1),1]);
mM = 2*eye(Num-1) - diag(ones(1,Num-2),1) - diag(ones(1,Num-2),-1);

% vx0 = [-1.502,-0.501,0.5,1.504];

vx0 = [-2,-0.4,0.4,2];
vx = vx0;
while 1
    
   tm = fn_m(vx);
   vFj = fn_Fj(vx);
   
   delta_vx = -tm\transpose(vFj);
   vx = vx+transpose(delta_vx);
   disp(vx);
   disp(norm(fn_Fj(vx)));
    
end

vx = vx0;

tm = fn_m(vx);

while 1
    vx_pre = vx;
    vL = zeros(1,Num);
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
        vL(p) = Num/2*f(x)-sum(f(x-tv));
        
    end
    
    disp(' vL = ');
    disp(vL);
    
    vN = vL(2:end) - vL(1:end-1);
    vN = vN';
    
    vD = -mM\vN;
    
    disp('vD = ');
    disp(vD);
    
    vd = 1./tanh(vD/2);
    
    disp('vd = ');
    disp(vd);
    
    x0 = 0 - sum(vd)/2;
    
    vx = [x0,vd']*m_d2x;
    
    disp('vx');
    disp(vx);
    
    if norm(vx_pre-vx)<1e-8
        break;
    end
    
end

while 1
    
    vx_total_pre = vx;
    vd = vx*m_x2d;
    for point = 1:(Num-1)
        while 1
            vx_pre = vx;
            x = vx(point);
            tv = vx;
            tv(point) = [];
            tv(point) = [];
            LHS = 4*Num/2*f(x);
            RHS = sum(f(x-tv));
            d = 1 + 2/(exp(LHS-RHS)-1);
            vd(point+1) = d; % d_point update !
            vx = vd*m_d2x;
            disp(point);
            disp(vx);
            disp(vd);
            disp(fn_Fj(vx));
            
            if norm(vx_pre-vx)<1e-8
                break;
            end
        end
        
    end
    
    
    vd_b = vx*m_x2d_b;
    
    for point = Num:-1:2
        while 1
            vx_pre = vx;
            x = vx(point);
            tv = vx;
            tv(point) = [];
            tv(point-1) = [];
            LHS = 4*Num/2*f(x);
            RHS = sum(f(x-tv));
            d = 1 + 2/(exp(-LHS+RHS)-1);
            vd_b(point-1) = d; % d_(point-1) update
            vx = vd_b*m_d2x_b;
            disp(point);
            disp(vx);
            disp(vd_b);
            disp(fn_Fj(vx));
            
            if norm(vx_pre-vx)<1e-8
                break;
            end
        end
        
    end
    
    disp(vx);
    disp(norm(vx_total_pre-vx));
    
    if norm(vx_total_pre-vx)<1e-8
        break;
    end
    
end


step = 0.1;
while step>1e-8
    vx_left = vx - step*ones(1,Num);
    vx_right = vx + step*ones(1,Num);
    
    value_left = sum(abs(fn_Fj(vx_left)));
    value_right = sum(abs(fn_Fj(vx_right)));
    value_center = sum(abs(fn_Fj(vx)));
    
    if value_center<min(value_left,value_right)
        step = step/2;
        disp(step);
    else
        if value_left<value_right
            vx = vx_left;
        else
            vx = vx_right;
        end
    end
    
    disp(vx);
    disp(sum(abs(fn_Fj(vx_right))));
    disp(vx(2:end)-vx(1:end-1));
    
end


function res = fn_Fj(vx)

Num = length(vx);
[mx1,mx2] = meshgrid(vx,vx);
tm = log(abs((mx1-mx2-1)./(mx1-mx2+1)));
tm(logical(eye(Num))) = 0;

vx1 = Num/2*log(abs((vx-1)./(vx+1)));
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