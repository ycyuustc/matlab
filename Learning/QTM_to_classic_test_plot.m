N = 8;
M = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%m_configuration = zeros(M,N);
m_configuration = floor(2*rand(M,N));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tv = 1:M;
tv = mod(tv,2)*2 - 1;
m_sign = diag(tv) * ones(M,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_op = m_configuration.*m_sign;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v_delta_x = [1,0,-1,0];
v_delta_y = [0,-1,0,1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_swerve = [2,1,4,3;1,2,3,4;4,3,2,1];



close(figure(1));
figure(1);
hold on;
axis([0,N,0,M]);
x_ori = 3; y_ori = 4; direction_ori = 1;
x = x_ori; y = y_ori; direction = direction_ori;
disp([x,y,direction]);
while 1  
   x_old = x;
   y_old = y;
   op = m_op(y,x); 
   direction = m_swerve(op+2,direction);
   
   dx = v_delta_x(direction);
   dy = v_delta_y(direction);
   
   x = x + dx;
   y = y + dy;
   
   x = mod(x-1,N) + 1;
   y = mod(y-1,M) + 1;
   
   disp([x,y,direction]);
   if abs(x_old - x)<N-1 && abs(y_old - y)<M-1
        plot([x_old,x],[y_old,y],'r*-');
   end
   
%    if x == x_ori && y == y_ori
%        break;
%    end
   
end

if direction_ori == direction
    disp("one_circle");
else
    disp("two_circle");
end





