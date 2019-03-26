close(figure(1));
figure(1);
hold on;
close(figure(2));
figure(2);
hold on;

a = 1;
b = 2;
[x,y] = fn_ana(a,b);

jude_distance = 0.35;

v_epsilon=[-7,-5,-3,-1,1,3,5,7];
v_state = zeros(1,length(v_epsilon));
g=0.01;

v_ind0=[2,3,5,6,8];

vx=v_epsilon(v_ind0);
vx=vx+g/2;
disp(vx);

v_ind = 1:length(v_ind0);

v_g = linspace(0.01,10,1000);

m_speed = zeros(length(vx),length(vx));
m_distance = fn_distance(vx);
m_distance_pre = m_distance;

for i=1:length(v_g)
   g = v_g(i);
   
   m_close_index = fn_judge(vx,judge_distance);
   m_distance_pre = m_distance;
   m_distance = fn_distance(vx);
   m_speed = m_distance-m_distance_pre;
   msp = m_speed;
   
   if isempty(m_close_index)
       vx = fn_steepest(vx,v_epsilon,v_ind,g);
   else
       [num_pair,~] = size(m_close_index);
       v_close_ind = reshape(m_close_index,1,numel(m_close_index));
       v_ind_t = v_ind;
       v_ind_t(v_close_ind) = [];
       
       while 1
           
           vx_pre = vx;
           for ipair = 1:num_pair
               ind2 = m_close_index(ipair,1);
               ind1 = m_close_index(ipair,2);
               x1 = vx(ind1); x2 = vx(ind2);
               x_center = (x1+x2)/2;
               [~,ind0] = min(abs(v_epsilon - x_center));
               
               vx = fn_steepest3(vx,v_epsilon,v_ind_t,g,ind1,ind2,ind0,msp);
               
           end
       
       if norm(vx_pre-vx)<1e-8
           break;
       end
       
       end
       
       
      
        
   end
   
  
   disp(vx);
   figure(1);
   plot(vx,g*ones(length(vx)),'r*');
   figure(2);
   plot(real(vx),imag(vx),'bo');
   pause(0.05);
   
   
   
end


vx2 = fn_solve2(vx,v_epsilon,g,ind1,ind2,ind0);
disp(vx2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%    Functions    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% solve 2 point, keep the rest unchanged. 
function res = fn_solve2(vx,v_epsilon,g,ind1,ind2,ind0)

vxr = vx; vxr([ind1,ind2]) = [];
vepr = v_epsilon; vepr(ind0) = [];
epsilon = v_epsilon(ind0);
x = vx(ind1) - epsilon; y = vx(ind2) - epsilon;
% lambda = y/x;
while 1
    x_pre = x; y_pre = y;
%     a = -fn_f(x+epsilon,vxr,vepr,g);
%     b = -fn_f(y+epsilon,vxr,vepr,g);
%     x = (a-b)/a/(a+b);
%     y = -(a-b)/b/(a+b);
    
    a = (x-y)*(fn_f(epsilon+x,vxr,vepr,g)-fn_f(epsilon+y,vxr,vepr,g))-2;
    lambda = a/2+sqrt(a^2-4)/2;
    [x,y] = fn_small_solve(lambda,vxr,vepr,g,epsilon);
    
    disp([x,y]);
    if norm([x,y]-[x_pre,y_pre])<1e-10
        break;
    end     
end

if real(x)+imag(x)<real(y)+imag(y)
    vx(ind1) = x + epsilon;
    vx(ind2) = y + epsilon;
else
    vx(ind1) = y + epsilon;
    vx(ind2) = x + epsilon;
end

res = vx;

end


function [x,y] = fn_small_solve(lambda,vxr,vepr,g,epsilon)

x = 0;
while 1
   x_pre = x;
   x = (1+lambda)/(1-lambda)/fn_f(epsilon+x,vxr,vepr,g);
   disp(x);
   if norm(x_pre - x)<1e-10
       break;
   end
end

y = lambda*x;

end

function res = fn_f(x,vxr,vepr,g)

res = sum(1./(x - vepr)) - 2*sum(1./(x - vxr)) - 2/g;

end

function v_result = fn_steepest3(vx0,v_epsilon,v_ind,g,ind1,ind2,ind0,msp)

vx = vx0;
while 1
    vx_pre = vx;
    vx = fn_steepest(vx,v_epsilon,v_ind,g);
    vx = fn_ana_2p(vx,v_epsilon,ind1,ind2,ind0,g,msp);
     
    if norm(vx_pre-vx)<1e-8
        break;
    end
    
end

v_result = vx;
end

function v_result = fn_steepest2(vx0,v_epsilon,v_ind,g,ind1,ind2,ind0)

vx = vx0;
while 1
    vx_pre = vx;
    vx = fn_steepest(vx,v_epsilon,v_ind,g);
    vx = fn_solve2(vx,v_epsilon,g,ind1,ind2,ind0);
     
    if norm(vx_pre-vx)<1e-8
        break;
    end
    
end

v_result = vx;
end

function v_result=fn_steepest(vx0,v_epsilon,v_ind,g)

% global diameter;
% global NeedAnalytic;
vx = vx0;
TotalStep=0;
while 1
    TotalStep=TotalStep+1;
    dFj=fn_Fj(vx,v_epsilon,g);
    Matrix=fn_m(vx,v_epsilon);
    dFj = dFj(v_ind);
    Matrix = Matrix(v_ind,v_ind);
    dvx=-((Matrix.')\(dFj.')).';
    
    vx(v_ind) = vx(v_ind) + dvx;
    
    tv = fn_Fj(vx,v_epsilon,g);
    
    TotalStep = TotalStep + 1;
    if norm(tv(v_ind))<1e-10
        break;
    end
      
end

v_result=vx;
disp(['totalstep=',num2str(TotalStep)]);

end


function v_result=fn_Fj(vx,v_epsilon,g)

M=length(vx);
[mj,mk]=meshgrid(vx,v_epsilon);
m_r=1./(mj-mk);
tv1=sum(m_r);

[mj,ml]=meshgrid(vx,vx);
m_r=(mj-ml);
m_r(logical(eye(M)))=ones(1,M);
m_r=1./m_r;
m_r(logical(eye(M)))=zeros(1,M);
tv2=sum(m_r);

v_result=tv1-2*tv2-2/g;

end


function m_result=fn_m(vx,v_epsilon)

M=length(vx);
[mj,mk]=meshgrid(vx,v_epsilon);
m_r=-1./(mj-mk).^2;
v_r=sum(m_r);
tm1=zeros(M,M);
tm1(logical(eye(M)))=v_r;

[mj,ml]=meshgrid(vx,vx);
m_r=(mj-ml).^2;
m_r(logical(eye(M)))=ones(1,M);
m_r=1./m_r;
m_r(logical(eye(M)))=zeros(1,M);
m_r=m_r*2;
v_r=sum(m_r);
tm2=zeros(M,M);
tm2(logical(eye(M)))=v_r;

tm3=-m_r;

m_result=tm1+tm2+tm3;

end



function vx = fn_ana_2p(vx,v_epsilon,ind1,ind2,ind0,g,msp)

[a,b] = fn_ab(vx,v_epsilon,ind1,ind2,ind0,g);
epsilon = v_epsilon(ind0);
x1 = vx(ind1); x2 = vx(ind2);
[x,y] = fn_ana(a,b);

if vs(ind0)==1
    if abs(real(x1-x2))>abs(imag(x1-x2))
        if abs(real(x-y))<abs(real(x1-x2))
            vx(ind1) = x + epsilon;
            vx(ind2) = y + epsilon;
        end
    else
        if abs(imag(x-y))>abs(real(x1-x2))
            vx(ind1) = x + epsilon;
            vx(ind2) = y + epsilon;
        end
    end
    
elseif vs(ind0) == -1
    if abs(real(x1-x2))>abs(imag(x1-x2))
        if abs(real(x-y))>abs(real(x1-x2))
            vx(ind1) = x + epsilon;
            vx(ind2) = y + epsilon;
        end
    else
        if abs(imag(x-y))<abs(real(x1-x2))
            vx(ind1) = x + epsilon;
            vx(ind2) = y + epsilon;
        end
    end  
    
end

% if abs(x-y)>0.45
%    if abs(real(x))>abs(imag(x))
%        x = -0.1;
%        y = 0.1;
%    else
%        x = real(x)-0.1j;
%        y = real(y)+0.1j;
%    end
% end
%     
% vx(ind1) = x + epsilon;
% vx(ind2) = y + epsilon;

end


function [x,y] = fn_ana(a,b)

c = sqrt(a^2+8*b);
if a<0
    delta = sqrt(2)*b*sqrt((a^2-4*b-a*c)/b^2);
    x = -(3*a+c-delta)/4/b;
    y = -(3*a+c+delta)/4/b;
else
    delta = sqrt(2)*b*sqrt((a^2-4*b+a*c)/b^2);
    x = -(3*a-c-delta)/4/b;
    y = -(3*a-c+delta)/4/b;
end

end

function [a,b] = fn_ab(vx,v_epsilon,ind1,ind2,ind0,g)

epsilon = v_epsilon(ind0);
vepr = v_epsilon;
vepr(ind0) = [];
vxr = vx;
vxr([ind1,ind2]) = [];

a = sum(1./(epsilon-vepr)) - 2*sum(1./(epsilon-vxr)) - 2/g;
b = -sum(1./(epsilon-vepr).^2)+ 2*sum(1./(epsilon-vxr).^2);

a = -a;
b = -b;

end


function res = fn_judge(vx,dis)

[mx1,mx2] = meshgrid(vx,vx);
dm = mx1 - mx2;

dm(triu(ones(length(vx)))) = 1;

[indx,indy] = find(dm<dis);

res = [indx,indy];

end

function dis = fn_distance(vx)

[mx1,mx2] = meshgrid(vx,vx);
dm = mx1-mx2;
dm(triu(ones(length(vx)))) = 0;
dis = dm;

end