% N*M size square lattice
M = 2; % y direction, or temperature direction
N = 3; % x direction  or site direction

% B = rand();
B = 0.1;
beta = 0.1; % minus temperature
lambda = beta/M;
b = 1/(1-lambda); % abs(b) on the permutation tensor
c = -lambda/(1-lambda);  % abs(c) on the identity tensor
logb = log(b);
logc = log(abs(c));
v_potential = [abs(c),abs(b)];

% res = fn_3_plus_2_exact(lambda + 0.0001,B) - fn_3_plus_2_exact(lambda,B);
% res = res/0.0001;
% disp(res);
[res1,res2,Z] = fn_3_plus_2_exact(lambda,B);
[res12,res22,Z2] = fn_3_plus_2_exact(lambda+0.0001,B);
fprintf('Z even = %f \n',res1);
fprintf('Z odd = %f \n',res2);
fprintf('Z odd over even = %f \n',res2/res1);
fprintf('cv frome Z even = %f \n',(res12-res1)/0.0001/res1);
fprintf('cv frome Z odd = %f \n',(res22-res2)/0.0001/res2);
fprintf('cv frome Z = %f \n',(Z2-Z)/0.0001/Z);

[res1,res2,Z] = fn_3_plus_2_exact(lambda,B);
[res12,res22,Z2] = fn_3_plus_2_exact(lambda,B+0.0001);
fprintf('B frome Z even = %f \n',(res12-res1)/0.0001/res1);
fprintf('B frome Z odd = %f \n',(res22-res2)/0.0001/res2);
fprintf('B frome Z = %f \n',(Z2-Z)/0.0001/Z);


m_sign = diag((mod(1:M,2)*2-1)) * ones(M,N);

% even configuration
m_configuration = ones(M,N);

% m_configuration = floor(rand(M,N)*2);
% m_configuration = zeros(M,N);
% m_configuration(1,1) = 1;

v_n = zeros(1,1+N);
v_n(1:N) = ones(1,N);
totCircle = N+M-1;

%% test code
% m_op = floor(rand(M,N)*3)-1;
% disp(m_op);
% 
% x = floor(rand()*N)+1;
% y = floor(rand()*M)+1;
% 
% disp([x,y]);
% m_op(y,x) = 1;
% [tn,tv] = circle_one_point(x,y,m_op);
% disp(tn);
% disp(tv);
% m_op(y,x) = 0;
% [tn,tv] = circle_one_point(x,y,m_op);
% disp(tn);
% disp(tv);
% m_op(y,x) = -1;
% [tn,tv] = circle_one_point(x,y,m_op);
% disp(tn);
% disp(tv);
% 
% B = 1;
% disp(fn_factor(x,y,m_op,B));
%%
Num_walk = 1000000;

for kkk = 1:100

m_v_n = zeros(Num_walk,N+1);
v_totCircle = zeros(1,Num_walk);
v_partial_lambda = zeros(1,Num_walk);
v_partial_B = zeros(1,Num_walk);
v_odd_over_even = zeros(1,Num_walk);
walk_step = 1;
while walk_step <= Num_walk
     
    v_n_pre = v_n;
    totCircle_pre = totCircle;

while 1
    x1 = floor(rand()*N)+1;
    y1 = floor(rand()*M)+1;
    x2 = floor(rand()*N)+1;
    y2 = floor(rand()*M)+1;
if x1~=x2 || y1~=y2
    break;
end
end

state1 = m_configuration(y1,x1);
state2 = m_configuration(y2,x2);
if state1 == 1 && state2 == 1
    
    flip_factor_1 = (c/b)^2;
elseif state1 == 0 && state2 == 0
    flip_factor_1 = (b/c)^2;
else
    flip_factor_1 = 1;
end

%% decision parameter 
m_configuration_1 = m_configuration;
m_configuration_1(y1,x1) = 1 - m_configuration_1(y1,x1);
[factor_0,num_circle_0,v0] = fn_factor(x1,y1,m_configuration.*m_sign,B);
[factor_1,num_circle_1,v1] = fn_factor(x1,y1,m_configuration_1.*m_sign,B);
factor_0_to_1 = factor_1/factor_0;
if num_circle_0 == 2 && num_circle_1 == 1
    % combine two circles
    pf = find(v_n==v0(1),1);
    v_n(pf) = v1(1);
    pf2 = find(v_n==v0(2),1);
    v_n(pf2) = 0;
    totCircle = totCircle - 1; 
end
if num_circle_0 == 1 && num_circle_1 == 2
   % seperate one circles
   pf = find(v_n==v0(1),1);
   v_n(pf) = v1(1);
   pf = find(v_n==0,1);
   v_n(pf) = v1(2);
   totCircle = totCircle + 1;
end

% disp('try step 1:');
% disp(v_n);

m_configuration_2 = m_configuration_1;
m_configuration_2(y2,x2) = 1 - m_configuration(y2,x2);
[factor_1,num_circle_1,v1] = fn_factor(x2,y2,m_configuration_1.*m_sign,B);
[factor_2,num_circle_2,v2] = fn_factor(x2,y2,m_configuration_2.*m_sign,B);
factor_1_to_2 = factor_2/factor_1;
if num_circle_1 == 2 && num_circle_2 == 1
    % combine two circles
    pf = find(v_n==v1(1),1);
    v_n(pf) = v2(1);
    pf2 = find(v_n==v1(2),1);
    v_n(pf2) = 0;
    totCircle = totCircle - 1; 
end
if num_circle_1 == 1 && num_circle_2 == 2
   % seperate one circles
   pf = find(v_n==v1(1),1);
   v_n(pf) = v2(1);
   pf = find(v_n==0,1);
   v_n(pf) = v2(2);
   totCircle = totCircle + 1;
end

% disp('try step 2:');
% disp(v_n);

flip_factor = flip_factor_1*factor_1_to_2*factor_0_to_1;
% disp(flip_factor);
%%
[factor_00,~,~] = fn_factor(1,1,m_configuration.*m_sign,B);
m_configuration_1 = m_configuration;
m_configuration_1(1,1) = 1 - m_configuration_1(1,1);
[factor_trans_11,~,~] = fn_factor(1,1,m_configuration_1.*m_sign,B);

index = m_configuration(1,1);
factor_11 = factor_trans_11/factor_00 ...
    *v_potential(1+(1-index))/v_potential(1+index);

if rand()>1/(1+flip_factor)
    % m_configuration 的值只在这里改变
    m_configuration(y1,x1) = 1 - m_configuration(y1,x1);
    m_configuration(y2,x2) = 1 - m_configuration(y2,x2);
else
    v_n = v_n_pre;
    totCircle = totCircle_pre;
end

partial_lambda = sum(sum(1/lambda/(1-lambda)-m_configuration./lambda));
v_partial_lambda(walk_step) = partial_lambda;
partial_B = sum(v_n.*tanh(v_n*B));
v_partial_B(walk_step) = partial_B;
m_v_n(walk_step,:) = v_n;
v_totCircle(walk_step) = totCircle;
v_odd_over_even(walk_step) = factor_11;

walk_step = walk_step + 1;

% if mod(walk_step,200000) == 0
%     disp(walk_step);
% end
% 


%% output data
% disp(flip_factor);
% disp(v_n);
% disp(totCircle);
% disp(m_configuration);
% disp('=========================');
%%


end

disp(kkk);
disp(mean(v_odd_over_even));

end

fprintf('cv from Monte Carlo = %f \n',mean(v_partial_lambda));
fprintf('M from Monte Carlo = %f \n', mean(v_partial_B));
fprintf('Odd over Even from Monte Carlo = %f \n', mean(v_odd_over_even));


% m_configuration_3 = m_configuration;
% m_configuration_3(y1,x1) = 1 - m_configuration_3(y1,x1);
% m_configuration_3(y2,x2) = 1 - m_configuration_3(y2,x2);
% 
% m_op = m_configuration_3.*m_sign;
% 
% c1 = ct{m_op(1,1)+2};
% c2 = ct{m_op(1,2)+2};
% c3 = ct{m_op(1,3)+2};
% c4 = ct{m_op(2,1)+2};
% c5 = ct{m_op(2,2)+2};
% c6 = ct{m_op(2,3)+2};
% c7 = expm([1,0;0,-1]*B);
% 
% res2 = ncon({c1,c2,c3,c4,c5,c6,c7,c7,c7},{[1,7,2,8],[2,10,3,11],[3,13,1,14],...
%     [4,8,5,9],[5,11,6,12],[6,14,4,15],[9,7],[12,10],[15,13]});
% 
% disp(res2/res1);
% disp(flip_factor);

function [res,num_circle,v_pearl] = fn_factor(x,y,m_op,B)
    [num_circle,v_pearl] = circle_one_point(x,y,m_op);
    if num_circle == 2
        res = 2*cosh(B*v_pearl(1))*2*cosh(B*v_pearl(2));
    else
        res = 2*cosh(B*v_pearl(1));
    end
end

function [num_circle,v_pearl] = circle_one_point(x,y,m_op)
[M,N] = size(m_op);
m_swerve = [2,1,4,3;1,2,3,4;4,3,2,1];
v_delta_x = [1,0,-1,0];
v_delta_y = [0,-1,0,1];

x_ori = x; y_ori = y;
direction_ori = 1;
op = m_op(y,x);

v_backward = [3,4,1,2];
direction_next = v_backward(m_swerve(op+2,direction_ori));

direction = 1;
which_circle = 1;
v_pearl = zeros(1,2);
break_flag = 0;
while 1
    if direction == 2 && y == M
        v_pearl(which_circle) = v_pearl(which_circle) + 1;
    end
    if direction == 4 && y == 1
        v_pearl(which_circle) = v_pearl(which_circle) + 1;
    end
    
    % update the leg, which is discribed by [x,y,direction]. 
    op = m_op(y,x);
    direction = m_swerve(op+2,direction); % update direction
    dx = v_delta_x(direction);
    dy = v_delta_y(direction);
    x = mod(x+dx-1,N)+1; % update x
    y = mod(y+dy-1,M)+1; % update y
      
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


function [res1,res2,Z] = fn_3_plus_2_exact(lambda,B)

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

b = 1/(1-lambda); % abs(b) on the permutation tensor
c = -lambda/(1-lambda);  % abs(c) on the identity tensor

c1 = b*m4p + c*m4id;
c4 = b*m4p2 + c*m4id;
c2 = c1; c3 = c1;
c5 = c4; c6 = c4;
c7 = expm([1,0;0,-1]*B);
Z = mcon({c1,c2,c3,c4,c5,c6,c7,c7,c7},...
    {[1,7,2,8],[2,10,3,11],[3,13,1,14],...
    [4,8,5,9],[5,11,6,12],[6,14,4,15],[9,7],[12,10],[15,13]});


c1 = cell(1,2);
c4 = cell(1,2);
c1{1} = b*m4p; c1{2} = c*m4id;
c2 = c1; c3 = c1;
c4{1} = b*m4p2; c4{2} = c*m4id;
c5 = c4; c6 = c4;
c7 = expm([1,0;0,-1]*B);

res1 = 0;
res2 = 0;
for i1 = 1:2
    for i2 = 1:2
        for i3 = 1:2
            for i4 = 1:2
                for i5 = 1:2
                    for i6 = 1:2
res = mcon({c1{i1},c2{i2},c3{i3},c4{i4},c5{i5},c6{i6},c7,c7,c7},...
    {[1,7,2,8],[2,10,3,11],[3,13,1,14],...
    [4,8,5,9],[5,11,6,12],[6,14,4,15],[9,7],[12,10],[15,13]});
if mod(i1+i2+i3+i4+i5+i6,2) == 0
    res1 = res1 + res;
else
    res2 = res2 + res;
end
                    end
                end
            end
        end
    end
end

end



