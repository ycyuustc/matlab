Num_walk = 10000;

beta = 0.1;
lambda = beta/M;
b = 1/(1-lambda);
c = -lambda/(1-lambda);
logb = log(b);
logc = log(abs(c));

M = 10;
N = 10;

% m_configuration = floor(2*rand(M,N));
m_sign = diag((mod(1:M,2)*2-1)) * ones(M,N);

tic;
circle_one_point(3,4,m_op);
toc;

m_configuration = ones(M,N);
m_op = m_configuration.*m_sign;
tic;
circle_one_point(3,4,m_op);
toc;

N_circle = N;
v_power_1 = zeros(1,Num_walk);
v_power_2 = zeros(1,Num_walk);
v_power_1(1) = M*N*logb;
v_power_2(1) = N*log(2);

for ir=2:Num_walk  
    x1 = floor(rand()*N) + 1;
    y1 = floor(rand()*M) + 1;
    while 1
        x2 = floor(rand()*N) + 1;
        y2 = floor(rand()*M) + 1;
        if x1~=x2 || y1~=y2
            break;
        end    
    end
    
    state1 = m_configuration(x1,y1);
    state2 = m_configuration(x2,y2);
    if state1 == 1 && state2 == 1
        flip_factor_1 = (c/b)^2;
    elseif state1 == 0 && state2 == 0
        flip_factor_1 = (b/c)^2;
    else
        flip_factor_1 = 1;
    end
    
    res1 = circle_one_point(x1,y1,m_op);
    if abs(m_op(x1,y1)) == 1
        m_op(x1,y1) = 0;
    else
        m_op(x1,y1) = m_sign(x1,y1);
    end
    res2 = circle_one_point(x2,y2,m_op);
    flip_factor_2 = res1 * res2;
    
    flip_factor = flip_factor_1 * flip_factor_2;
    
    if rand()<flip_factor/(1+flip_factor)
       m_configuration(x1,y1) = 1-m_configuration(x1,y1);
       m_configuration(x2,y2) = 1-m_configuration(x2,y2);
       v_power_1(ir) = v_power_1(ir-1) + log(flip_factor_1);
       v_power_2(ir) = v_power_2(ir-1) + log(flip_factor_2);
    else
        v_power_1(ir) = v_power_1(ir-1);
        v_power_2(ir) = v_power_2(ir-1);
    end
    m_op = m_configuration.*m_sign;
    disp(ir);
    disp(flip_factor_1);
    disp(flip_factor_2);

end

function res = circle_one_point(x,y,m_op)
[M,N] = size(m_op);
m_swerve = [2,1,4,3;1,2,3,4;4,3,2,1];
v_delta_x = [1,0,-1,0];
v_delta_y = [0,-1,0,1];

x_ori = x; y_ori = y;
direction_ori = 1;
direction = 1;


while 1  
   op = m_op(y,x); 
   direction = m_swerve(op+2,direction);  
   dx = v_delta_x(direction);
   dy = v_delta_y(direction);
   x = mod(x+dx-1,N)+1;
   y = mod(y+dy-1,M)+1;  
  
   if x == x_ori && y == y_ori
       break;
   end
   
end

if direction_ori == direction
    res = 1/2;
else
    res = 2;
end

end







