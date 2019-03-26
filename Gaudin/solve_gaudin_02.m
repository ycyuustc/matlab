close(figure(1));
figure(1);
hold on;
close(figure(2));
figure(2);
hold on;

a = 1;
b = 2;
[x,y] = fn_ana(a,b);

dis = 0.25;

v_epsilon=[-7,-5,-3,-1,1,3,5,7];
v_state = zeros(1,length(v_epsilon));
g=0.01;

v_ind0=[2,3,5,6,8];

vx=v_epsilon(v_ind0);
vx=vx+g/2;
vx_pre = vx;
disp(vx);

v_ind = 1:length(v_ind0);
v_epsilon_state = zeros(1,length(vx));

v_g = linspace(0.01,10,1000);
g = 0.01;
g_step = 0.01;
while 1
   % update the state
   v_epsilon_state = fn_state(vx,vx_pre,v_epsilon,dis,v_epsilon_state);
   [m_close_index,v_ind0] = fn_find_index(vx,v_epsilon,v_epsilon_state,dis);
   vx_pre = vx;
   if isempty(v_ind0)
       vx = fn_steepest(vx,v_epsilon,v_ind,g);
       g = g+g_step;
   else
       
       v_ind_t = v_ind;
       v_close_index = reshape(m_close_index,1,numel(m_close_index));
       v_ind_t(v_close_index) = [];
       
%        while 1
           
           vx_pre2= vx;
           for ipair = 1:length(v_ind0)
               ind2 = m_close_index(ipair,1);
               ind1 = m_close_index(ipair,2);
               ind0 = v_ind0(ipair);       
               vx = fn_steepest3(vx,v_epsilon,v_ind_t,g,...
                   ind1,ind2,ind0,v_epsilon_state);
               disp(vx);
           end
       
%        if norm(vx_pre2-vx)<1e-8
%            break;
%        end
       
%        end
       
       g = g + g_step;
      
        
   end
   
  
   disp(vx);
   figure(1);
   plot(vx,g*ones(length(vx)),'r*');
   figure(2);
   plot(real(vx),imag(vx),'bo');
   pause(0.05);
   
   if g > 20
       break;
   end
   
   
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

function v_result = fn_steepest3(vx0,v_epsilon,v_ind,g,ind1,ind2,ind0,v_epsilon_state)

vx = vx0;
% while 1
%     vx_pre = vx;
%     vx = fn_steepest(vx,v_epsilon,v_ind,g);
    vx = fn_ana_2p(vx,v_epsilon,ind1,ind2,ind0,g,v_epsilon_state);
%     vx = fn_steepest(vx,v_epsilon,v_ind,g);
%     if norm(vx_pre-vx)<1e-8
%         break;
%     end
%     
% end

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



function vx = fn_ana_2p(vx,v_epsilon,ind1,ind2,ind0,g,v_epsilon_state)

[a,b] = fn_ab(vx,v_epsilon,ind1,ind2,ind0,g);
epsilon = v_epsilon(ind0);

[x,y] = fn_ana(a,b);
z0 = v_epsilon_state(ind0);
z = x - y;

if abs(real(z0*conj(z)))<abs(z0)^2
   x1 = x + epsilon;
   x2 = y + epsilon;
   
   if abs(vx(ind1)-x1)+abs(vx(ind2)-x2)<abs(vx(ind1)-x2)+abs(vx(ind2)-x1)
       vx(ind1) = x1;
       vx(ind2) = x2;
   else
       vx(ind1) = x2;
       vx(ind2) = x1;
   end

end



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

function v_epsilon_state = fn_state(vx,vx_pre,v_epsilon,dis,v_epsilon_state)


num = length(v_epsilon);
v_act_pre = fn_activate(vx_pre,v_epsilon,dis);
v_act = fn_activate(vx,v_epsilon,dis);
for i = 1:num
   if v_act(i)==1 && v_act_pre(i)==0
       epsilon = v_epsilon(i);
       tv = vx - epsilon;
       v_ind = find(abs(tv)<dis);
       ind1 = v_ind(1);
       ind2 = v_ind(2);
       x1 = vx(ind1);
       x2 = vx(ind2);
       v_epsilon_state(i) = x1 - x2;
   end
   
end

v_act_pre = fn_activate(vx_pre,v_epsilon,dis*1.2);
v_act = fn_activate(vx,v_epsilon,dis*1.2);
for i = 1:num
    if v_act(i)==0 && v_act_pre(i)==1
       v_epsilon_state(i) = 0;
   end
end

end

function v_activate = fn_activate(vx,v_epsilon,dis)

num = length(v_epsilon);
v_activate = zeros(1,num);

for i = 1:num
   
    epsilon = v_epsilon(i);
    tv = vx - epsilon;
    v_ind = find(abs(tv)<dis);
    if length(v_ind)>=2
        v_activate(i) = 1;
    end
    
end

end

function [m_close_ind,v_ind0] = fn_find_index(vx,v_epsilon,v_epsilon_state,dis)

v_ind0 = find(v_epsilon_state~=0);
num = length(v_ind0);
m_close_ind = zeros(num,2);

for i=1:num
    
    ind0 = v_ind0(i);
    epsilon = v_epsilon(ind0);
    tv = find(abs(vx - epsilon)<dis*1.2);
    m_close_ind(i,:) = tv;
    
end

end