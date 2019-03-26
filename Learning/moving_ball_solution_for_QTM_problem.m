withfigure=0;

if withfigure==1
    close(figure(1));
    figure(1);
    hold on;
    close(figure(2));
    figure(2);
    hold on;
    close(figure(3));
    figure(3);
    hold on;
end


% vx = zeros(1,Num);
% 
vx = zeros(1,8);
vx(1) = -6;
vx(2) = -4;
vx(3) = -2;
vx(4) = -0.51;
vx(5) = 0.51;
vx(6) = 2;
vx(7) = 4;
vx(8) = 6;
vx(7) = 2.53;
vx(8) = 3.54;
% tv1 = -10.5:10.5;
% Num = length(tv1);
% vx0 = linspace(-0.1,0.1,Num)+tv1;
Num = length(vx);

global m_connection;
m_connection = eye(Num);


% vx = vx0;

vv = zeros(1,Num);

v_force = zeros(1,Num);
y = 0;

count = 0;
while 1
    
    for i=1:Num
        v_force(i) = fn_force(vx,i);
    end
    
    delta_T = 1e-4;
    mass = 4.0;
    va = v_force/mass;
    mu2 = 0.05;
    va = va  - 0.2*(exp(abs(vv)/mu2)-1).*sign(vv);
    
    
    vx_pre = vx;
    vx = vx + delta_T*vv;
    vv = vv + delta_T*va;
    
    [vx,vv,va] = fn_adjust(vx,vv,va,vx_pre);
    
    y = y+delta_T;
    
    disp(vx);
    disp(vv);
    disp(va);
    
    if withfigure==1
        figure(1);
        plot(vx(1),y,'b*');
        plot(vx(2),y,'g*');
        plot(vx(3),y,'r*');
        plot(vx(4),y,'y*');
        %     plot(vx(5),y,'b*');
        %     plot(vx(6),y,'g*');
        figure(2);
        plot(vv(1),y,'b*');
        plot(vv(2),y,'g*');
        plot(vv(3),y,'r*');
        plot(vv(4),y,'y*');
        %     plot(vv(5),y,'b*');
        %     plot(vv(6),y,'g*');
        figure(3);
        plot(va(1),y,'b*');
        plot(va(2),y,'g*');
        plot(va(3),y,'r*');
        plot(va(4),y,'y*');
        %     plot(va(5),y,'b*');
        %     plot(va(6),y,'g*');
        
        pause(0.01);
        
    end
    
    if norm(va)<1e-6
        break;
    end
    
    count = count + 1;
    if mod(count,1000) == 0
        disp('have a rest');
    end
    
end

function res = fn_force(vx,ind)
global m_connection
v_judge = m_connection(ind,:);
mass = sum(v_judge);
v_connect = find(v_judge==1);
temp_sum = 0;
for i=1:length(v_connect)
    temp_sum = temp_sum + fn_force_onepoint(vx,v_connect(i));
end
res = temp_sum/mass;
end

function res = fn_force_onepoint(vx,ind)
global m_connection;
f = @(x) log(abs(x-1)./abs(x+1));
Num = length(vx);

v_judge = m_connection(ind,:);
tv = vx(logical(v_judge));
mass = sum(v_judge);
LHS = 2*Num/2*sum(f(tv));

v_judge = 1-v_judge;
tv = vx(logical(v_judge));
RHS = sum(f(vx(ind)-tv));

res = (LHS - RHS)/mass;

end

function [vx,vv,va] = fn_adjust(vx,vv,va,vx_pre)
global m_connection;
Num = length(vx);

m1 = ones(Num,1)*vx;
m2 = transpose(m1);
m_sign_1 = sign(abs(m1-m2)-1);

m1 = ones(Num,1)*vx_pre;
m2 = transpose(m1);
m_sign_2 = sign(abs(m1-m2)-1);

[ind_1,ind_2] = find(m_sign_1~=m_sign_2);

for i = 1:length(ind_1)
   if m_connection(ind_1(i),ind_2(i))==0
   m_connection(ind_1(i),ind_2(i))=1; 
   p1 = ind_1(i);
   p2 = ind_2(i);
   x1 = vx(p1);
   x2 = vx(p2);
   if x1<x2
      v_ind_left = find(vx<=x1);
      v_ind_right = find(vx>=x2);
      move_distance = 1-(x2-x1);
      vx(v_ind_left) = vx(v_ind_left) - move_distance/2;
      vx(v_ind_right) = vx(v_ind_right) + move_distance/2;
   else
      v_ind_left = find(vx<=x2);
      v_ind_right = find(vx>=x1);
      move_distance = 1-(x1-x2);
      vx(v_ind_left) = vx(v_ind_left) - move_distance/2;
      vx(v_ind_right) = vx(v_ind_right) + move_distance/2;      
   end
%    v_mean = (vv(p1)+vv(p2))/2;
%    a_mean = (va(p1)+va(p2))/2;
%    vv(p1) = v_mean;
%    vv(p2) = v_mean;
%    va(p1) = a_mean;
%    va(p2) = a_mean;
   end
end

while 1
    m_connect_pre = m_connection;
    m_connection = sign(m_connection*m_connection);
    if isequal(m_connect_pre, m_connection)
        break;
    end
end

% for i=1:Num-1
%     if vx(i+1)-vx(i)<=1
%         left_move_distance = 1-(vx(i+1)-vx(i));
%         vx(1:i) = vx(1:i) - left_move_distance;
%         v_mean = (vv(i+1)+vv(i))/2;
%         vv(i+1) = v_mean;
%         vv(i) = v_mean;
%         a_mean = (va(i+1)+va(i))/2;
%         va(i+1) = a_mean;
%         va(i) = a_mean;
%         m_connection(i,i+1) = 1;
%         m_connection(i+1,i) = 1;
%         while 1
%             m_connect_pre = m_connection;
%             m_connection = sign(m_connection*m_connection);
%             if isequal(m_connect_pre, m_connection)
%                 break;
%             end
%         end
%     end
% end

for i = 1:Num
   v_judge = m_connection(i,:);
   v_connect = find(v_judge==1);
   mass = sum(v_judge);
   sum_v = 0;
   sum_a = 0;
   for k = 1:length(v_connect)
       sum_v = sum_v + vv(v_connect(k));
       sum_a = sum_a + va(v_connect(k));
   end
   for k = 1:length(v_connect)
       vv(v_connect(k)) = sum_v/mass;
       va(v_connect(k)) = sum_a/mass;
   end
end

end

