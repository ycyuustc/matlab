close(figure(1));
figure(1);
hold on;
close(figure(2));
figure(2);
hold on;

% fn_onepoint(vx,v_epsilon,g,index);


v_epsilon=[1,3,5,7,9,11];
v_state = zeros(1,length(v_epsilon));
g=0.01;
b = 0.1;
dis = 0.2;

v_ind0=[1,2,3,5];

vx=v_epsilon(v_ind0);
vx=vx+g/2;
v_ind = 1:length(vx);

vx1 = fn_onepoint(vx,v_epsilon,g,b,2);
disp(fn_Fj(vx1,v_epsilon,g,b));

g_step = 0.01;
while 1
    [m_close_index,v_close_index] = fn_judge(vx,v_epsilon,dis);
    
    if isempty(v_close_index)
        vx = fn_steepest(vx,v_epsilon,g,v_ind);  
    else
        v_ind_t = v_ind;
        tv = reshape(m_close_index,1,numel(m_close_index));
        v_ind_t(tv) = [];
        vx = fn_steepest2(vx,v_epsilon,g,v_ind_t,m_close_index,v_close_index);
    end
    disp(vx);
    disp(norm(fn_Fj(vx,v_epsilon,g)));
    disp(g);
    figure(1);
    plot(vx,g*ones(length(vx)),'r*');
    figure(2);
    plot(real(vx),imag(vx),'bo');
    pause(0.05);
    
    
    g = g + g_step;
end

disp(vx);
disp(fn_Fj(vx,v_epsilon,g));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [m_close_indexa,v_inda] = fn_judge(vx,v_epsilon,dis)

[mx1,mx2] = meshgrid(vx,vx);
dm = mx1-mx2;
dm(logical(triu(ones(length(vx))))) = 1;
dm = abs(dm);
[vinda2,vinda1] = find(dm<dis);

[mx1,mx2] = meshgrid(vx,vx);
dm = mx1 + mx2;
dm(logical(triu(ones(length(vx))))) = 1;
dm = abs(dm);
[vindb2,vindb1] = find(dm<dis);


m_close_indexa = [vinda1,vinda2];
v_inda = zeros(1,length(vinda1));
for i = 1:length(vind1)
   ind1 = vind1(i);
   ind2 = vind2(i);
   x1 = vx(ind1);
   x2 = vx(ind2);
   x_center = (x1+x2)/2;
   tv = v_epsilon - x_center;
   [~,ind0] = min(abs(tv));
   v_inda(i) = ind0;
end

end


function vx = fn_steepest2(vx,v_epsilon,g,v_ind_t,m_close_index,v_close_index)

while 1
   vx_pre = vx;
   for k=1:length(v_close_index)
      ind1 = m_close_index(k,1);
      ind2 = m_close_index(k,2);
      ind0 = v_close_index(k);
      vx = fn_ana_2p(vx,v_epsilon,ind1,ind2,ind0,g);
   end
%    vx = fn_steepest3(vx,v_epsilon,g,v_ind_t);
   if norm(vx_pre-vx)<1e-8
       break;
   end
end

end

function vx = fn_steepest(vx,v_epsilon,g,v_ind)

while 1
   vx_pre = vx;
   for i = 1:length(v_ind)
       vx = fn_onepoint(vx,v_epsilon,g,v_ind(i));
   end
    
   disp(vx);
   if norm(vx_pre - vx)<1e-6
       break;
   end
    
end


end


function vx = fn_steepest3(vx,v_epsilon,g,v_ind)

while 1
   vx_pre = vx;
   for i = 1:length(v_ind)
       vx = fn_onepoint(vx,v_epsilon,g,v_ind(i));
   end
    
   disp(vx);
   if norm(vx_pre - vx)<1e-6
       break;
   end
    
end


end


function vx = fn_onepoint(vx,v_epsilon,g,b,index)

vxr = vx;
vxr(index) = [];

x = vx(index);
while 1
    
   f = sum(1./(x-v_epsilon))+sum(1./(x+v_epsilon))...
       - 2*sum(1./(x-vxr)) -2*sum(1./(x+vxr))...
       - 2/g*(1/(x-1j*b)+1/(x+1j*b));
   df = -sum(1./(x-v_epsilon).^2) -sum(1./(x-v_epsilon).^2)...
        + 2*sum(1./(x-vxr).^2)+2*sum(1./(x+vxr).^2)...
        +2/g*(1/(x-1j*b)^2+1/(x+1j*b)^2);
   
   x = x - f/df;
   
   if abs(f)<1e-10
       break;
   end
    
end

vx(index) = x;

end


function v_result=fn_Fj(vx,v_epsilon,g,b)

M=length(vx);
[mj,mk]=meshgrid(vx,v_epsilon);
m_r=1./(mj-mk)+1./(mj+mk);
tv1=sum(m_r);

[mj,ml]=meshgrid(vx,vx);
m_r=(mj-ml);
m_r(logical(eye(M)))=ones(1,M);
m_r=1./m_r;
m_r(logical(eye(M)))=zeros(1,M);
tv2=sum(m_r);

m_r3 = (mj+ml);
m_r3(logical(eye(M)))=ones(1,M);
m_r3=1./m_r3;
m_r3(logical(eye(M)))=zeros(1,M);
tv3=sum(m_r3);


v_result=tv1-2*tv2-2*tv3-2/g*(1./(vx-1j*b)+1./(vx+1j*b));

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



function vx = fn_ana_2p(vx,v_epsilon,ind1,ind2,ind0,g)

[a,b] = fn_ab(vx,v_epsilon,ind1,ind2,ind0,g);
epsilon = v_epsilon(ind0);

[x,y] = fn_ana(a,b);

x1 = x + epsilon;
x2 = y + epsilon;

if abs(vx(ind1)-x1)+abs(vx(ind2)-x2)<abs(vx(ind1)-x2)+abs(vx(ind2)-x1)
    vx(ind1) = x1;
    vx(ind2) = x2;
else
    vx(ind1) = x2;
    vx(ind2) = x1;
end


% z0 = v_epsilon_state(ind0);
% z = x - y;
% 
% if abs(real(z0*conj(z)))<abs(z0)^2
%    x1 = x + epsilon;
%    x2 = y + epsilon;
%    
%    if abs(vx(ind1)-x1)+abs(vx(ind2)-x2)<abs(vx(ind1)-x2)+abs(vx(ind2)-x1)
%        vx(ind1) = x1;
%        vx(ind2) = x2;
%    else
%        vx(ind1) = x2;
%        vx(ind2) = x1;
%    end
% 
% end

end