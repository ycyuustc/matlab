M = 40;
N = 4;

B = 0.0;
beta = 0.1;
lambda = beta/M;
theta = lambda/(1-lambda);

%%%%%%%%% test code %%%%%%%%%%%%%%%%%%%%%%%%%%
m_configuration = sign(rand(M,N)-0.5);
m_op = m_configuration;
x = floor(rand()*N) + 1;
y = floor(rand()*M) + 1;
[res,num_circle,m_pearl] = fn_factor(x,y,m_op,B);
disp(m_op);
fprintf('x = %d ; y = %d \n',x,y);
fprintf('res = %f\n',res);
fprintf('number of circles = %d\n',num_circle);
fprintf('pearl matrix: \n');
disp(m_pearl);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z1 = fn_original(lambda,B,M,N);
z2 = fn_original(lambda+1e-6,B,M,N);
partial_theta_exact = (1-lambda)^2*(z2-z1)/z1/1e-6;

z1 = fn_original(lambda,B,M,N);
z2 = fn_original(lambda,B+1e-6,M,N);
partial_B_exact = (z2-z1)/z1/1e-6;



Num_poach = 100;
Num_walk = 100000;

close(figure(1)); figure(1); hold on;
plot(1:Num_poach,partial_theta_exact*ones(1,Num_poach));

% initialization 
m_configuration = -diag((mod(1:M,2)*2-1)) * ones(M,N);
m_n = zeros(2,N);
m_n(1,:) = mod(1:N,2);
m_n(2,:) = mod(2:(N+1),2);
totCircle = N;
v_partial_theta = zeros(1,Num_walk);
v_partial_B = zeros(1,Num_walk);

v_theta = zeros(1,Num_poach);
v_B = zeros(1,Num_poach);

for ipoach = 1:Num_poach
for kkk = 1:Num_walk
    m_n_pre = m_n;
    totCircle_pre = totCircle;
    
    x1 = floor(rand()*N)+1;
    y1 = floor(rand()*M)+1;
    
    m_configuration_1 = m_configuration;
    m_configuration_1(y1,x1) = - m_configuration(y1,x1);
    [factor_0,num_circle_0,m_0] = fn_factor(x1,y1,m_configuration,B);
    [factor_1,num_circle_1,m_1] = fn_factor(x1,y1,m_configuration_1,B);
    factor_0_to_1 = factor_1/factor_0;
    
    spin = m_configuration(y1,x1);
    local_0 = (1+theta)/2 + (1-theta)/2*spin*(-1)^y1;
    local_1 = (1+theta)/2 + (1-theta)/2*(-spin)*(-1)^y1;
    local_0_to_1 = local_1/local_0;
    
    flip_factor = factor_0_to_1*local_0_to_1;
    
    % check: the flip factor by exactly trace the tensor network
%     factor_exact_0 = fn_4_plus_2_exact2(theta,B,m_configuration);
%     factor_exact_1 = fn_4_plus_2_exact2(theta,B,m_configuration_1);
  
%     flip_factor_exact = factor_exact_1/factor_exact_0;
%     disp(flip_factor);
%     disp(flip_factor_exact);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if rand()>1/(1+flip_factor)
        m_configuration(y1,x1) = - m_configuration(y1,x1);
        if num_circle_0 == 2 && num_circle_1 == 1
            % combine two circles
            pf = find(m_n(1,:)==m_0(1,1)&m_n(2,:)==m_0(2,1),1);
            m_n(:,pf) = m_1(:,1);
            pf2 = find(m_n(1,:)==m_0(1,2)&m_n(2,:)==m_0(2,2),1);
            m_n(:,pf2) = [0;0];
            totCircle = totCircle - 1;
        end
        if num_circle_0 == 1 && num_circle_1 == 2
            % seperate one circles
            pf = find(m_n(1,:)==m_0(1,1)&m_n(2,:)==m_0(2,1),1);
            m_n(:,pf) = m_1(:,1);
            pf = find(m_n(1,:)==0&m_n(2,:)==0,1);
            m_n(:,pf) = m_1(:,2);
            totCircle = totCircle + 1;
        end
        
    end
    
%     % display
%     fprintf('configuration matrix: \n');
%     disp(m_configuration);
%     fprintf('pearl matrix: \n');
%     disp(m_n);
    % end display
    partial_theta = sum(sum(1-diag(power(-1,1:M))*m_configuration))...
        /2/theta;
    v_partial_theta(kkk) = partial_theta;
    v_rb = m_n(1,:) - m_n(2,:);
    partial_B = sum(v_rb.*tanh(v_rb*B));
    v_partial_B(kkk) = partial_B;
    
%     if mod(kkk,10000) == 0
%         disp(kkk);
%     end
    
end

fprintf('the exact partial theta = %f \n',partial_theta_exact);
fprintf('the MC partial theta = %f \n',mean(v_partial_theta(1:end)));

fprintf('the exact partial B = %f \n',partial_B_exact);
fprintf('the MC partial B = %f \n',mean(v_partial_B(1:end)));

v_theta(ipoach) = mean(v_partial_theta);
v_B(ipoach) = mean(v_partial_B);

figure(1);
hold on;
plot(ipoach,v_theta(ipoach),'ro');
plot(ipoach,mean(v_theta(1:ipoach)),'b*');
pause(0.1);


end




disp('End');



function [res,num_circle,m_pearl] = fn_factor(x,y,m_op,B)
[num_circle,m_pearl] = circle_one_point(x,y,m_op);
if num_circle == 2
    res = 2*cosh(B*(m_pearl(1,1)-m_pearl(2,1)))...
        *2*cosh(B*(m_pearl(1,2)-m_pearl(2,2)));
else
    res = 2*cosh(B*(m_pearl(1,1)-m_pearl(2,1)));
end

end


function [num_circle,m_pearl] = circle_one_point(x,y,m_op)
[M,N] = size(m_op);
m_swerve = [4,3,2,1;1,2,3,4;2,1,4,3];
v_delta_x = [1,0,-1,0];
v_delta_y = [0,-1,0,1];

x_ori = x; y_ori = y;
direction_ori = 1;
op = m_op(y,x);

v_backward = [3,4,1,2];
direction_next = v_backward(m_swerve(op+2,direction_ori));

direction = 1;
which_circle = 1;
m_pearl = zeros(2,2);
break_flag = 0;
while 1
    if (direction == 2 && y == M) || (direction == 4 && y == 1)
        if mod(x,2) ==1
            % red pearl
            m_pearl(1,which_circle) = m_pearl(1,which_circle) + 1;
        else
            % blue pearl
            m_pearl(2,which_circle) = m_pearl(2,which_circle) + 1;
        end
    end
    
    % update the leg, which is discribed by [x,y,direction].
    op = m_op(y,x);
    direction = m_swerve(op+2,direction); % update direction
    dx = v_delta_x(direction);
    dy = v_delta_y(direction);
    x = mod(x+dx-1,N)+1; % update x
    y = mod(y+dy-1,M)+1; % update y
    %     disp([x,y,direction]);
    
    %  finished the update. now [x,y,direction] is the new leg.
    if x == x_ori && y == y_ori  % if go back to the original point
        if direction == direction_ori
            if break_flag == 1
                break;
            else
                direction_ori = ...
                    find((1:4)~=direction_ori&(1:4)~=direction_next,1);
                direction = direction_ori;
                which_circle = 2;
                break_flag = 1;
            end
        else
            break_flag = 1;
        end
    end
    
end  % end circle running
num_circle = which_circle;

end




function Z = fn_4_plus_2_exact(theta,B)

m4id=zeros(2,2,2,2);
m4id(1,1,1,1)=1;
m4id(1,2,1,2)=1;
m4id(2,1,2,1)=1;
m4id(2,2,2,2)=1;

m4p=zeros(2,2,2,2);
m4p(1,1,1,1)=1;
m4p(1,2,2,1)=1;
m4p(2,1,1,2)=1;
m4p(2,2,2,2)=1;

m4p2 = permute(m4p,[1,4,3,2]);

R1 = m4p + theta*m4p2;
R2 = m4p2 + theta*m4p;

mB1 = expm([1,0;0,-1]*B);
mB2 = expm(-[1,0;0,-1]*B);
Z = mcon({R1,R1,R1,R1,R2,R2,R2,R2,mB1,mB2,mB1,mB2},...
    {[1,10,2,11],[2,13,3,14],[3,16,4,17],[4,19,1,20],...
    [5,11,6,9],[6,14,7,12],[7,17,8,15],[8,20,5,18],...
    [9,10],[12,13],[15,16],[18,19]});

end


function Z = fn_4_plus_2_exact2(theta,B,m_configuration)

m4id=zeros(2,2,2,2);
m4id(1,1,1,1)=1;
m4id(1,2,1,2)=1;
m4id(2,1,2,1)=1;
m4id(2,2,2,2)=1;

m4p=zeros(2,2,2,2);
m4p(1,1,1,1)=1;
m4p(1,2,2,1)=1;
m4p(2,1,1,2)=1;
m4p(2,2,2,2)=1;

m4p2 = permute(m4p,[1,4,3,2]);

m_op = m_configuration;

% R1 = m4p + theta*m4p2;
% R2 = m4p2 + theta*m4p;

mB1 = expm([1,0;0,-1]*B);
mB2 = expm(-[1,0;0,-1]*B);

tc = cell(3,4);
for y = 1:2
    for x = 1:4
        tc{y,x} = (1+m_op(y,x))/2*m4p2+(1-m_op(y,x))/2*m4p;
        if (-1)^y*m_op(y,x)<0
            tc{y,x} = tc{y,x}*theta;
        end
    end
end

for x = 1:4
   if mod(x,2) == 1
       tc{3,x} = mB1;
   else
       tc{3,x} = mB2;
   end
end

tc = reshape(transpose(tc),[1,3*4]);
Z = mcon(tc,...
    {[1,10,2,11],[2,13,3,14],[3,16,4,17],[4,19,1,20],...
    [5,11,6,9],[6,14,7,12],[7,17,8,15],[8,20,5,18],...
    [9,10],[12,13],[15,16],[18,19]});

end


function Z = fn_original(lambda,B,M,N)

m4id=zeros(2,2,2,2);
m4id(1,1,1,1)=1;
m4id(1,2,1,2)=1;
m4id(2,1,2,1)=1;
m4id(2,2,2,2)=1;

m4p=zeros(2,2,2,2);
m4p(1,1,1,1)=1;
m4p(1,2,2,1)=1;
m4p(2,1,1,2)=1;
m4p(2,2,2,2)=1;

m4p2 = permute(m4p,[1,4,3,2]);

R1 = 1/(1-lambda)*m4p + lambda/(lambda-1)*m4id;
R2 = 1/(1-lambda)*m4p2 + lambda/(lambda-1)*m4id;

tct = cell(1,N);
tci = cell(1,N);
for k = 1:(N-1)
   tct{k} = R1;
   tci{k} = [k,-k,k+1,-k-N];
end
tct{N} = R1;
tci{N} = [N,-N,1,-N-N];
T1 = mcon(tct,tci);

tct = cell(1,N);
tci = cell(1,N);
for k = 1:(N-1)
   tct{k} = R2;
   tci{k} = [k,-k,k+1,-k-N];
end
tct{N} = R2;
tci{N} = [N,-N,1,-N-N];
T2 = mcon(tct,tci);

T12 = mcon({T1,T2},{[-(1:N),(1:N)],[(1:N),-(N+(1:N))]});

tensor = T12;
for k = 1:(M/2-1)
   tensor = mcon({tensor,T12},{[-(1:N),(1:N)],[(1:N),-(N+(1:N))]}); 
end

mB = expm(B*[1,0;0,-1]);
tct = cell(1,N);
tci = cell(1,N);
for k = 1:N
   tct{k} = mB;
   tci{k} = [-k,-(N+k)];
end
tensor_B = ncon(tct,tci);

Z = mcon({tensor_B,tensor},{[(1:N),(N+(1:N))],[(N+(1:N)),(1:N)]});

end