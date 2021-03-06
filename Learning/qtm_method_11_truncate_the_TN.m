global max_error; max_error = 0;
M = 1600;
N = 10;

B = 0.0;
beta = 1.0;
lambda = beta/M;
theta = lambda/(1-lambda);

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

tensornet = cell(M,N);
for m=1:M
    for n=1:N
        if mod(m,2)==1
            tensornet{m,n}=R1;
        else
            tensornet{m,n}=R2;
        end
    end
end


tv1 = 1:N;
tv2 = -tv1;
tv3 = [2:N,1];
tv4 = -tv1-N;
tm = [tv1;tv2;tv3;tv4];
tm = transpose(tm);
cind = cell(1,N);
for k=1:N
   cind{k} = tm(k,:); 
end
ct1 = tensornet(1,:);
ct2 = tensornet(2,:);
m1 = mcon(ct1,cind);
m2 = mcon(ct2,cind);
m1 = reshape(m1,2^N,2^N);
m2 = reshape(m2,2^N,2^N);
% [cTr,v_factor] = fn_add_line(tensornet(1,:),tensornet(2,:),D);

% test function fn_add_line
t = cell(2,4);
for x=1:2
    for y=1:4
        t{x,y} = rand(2,2,2,2);
    end
end
% t = tensornet;
% [tr,tf] = fn_add_line(t(1,:),t(2,:),4);
% Z0 = mcon({t{1,1},t{1,2},t{1,3},t{1,4},t{2,1},t{2,2},t{2,3},t{2,4}},...
%     {[1,9,2,10],[2,11,3,12],[3,13,4,14],[4,15,1,16],[5,10,6,9],...
%     [6,12,7,11],[7,14,8,13],[8,16,5,15]});
% Z1 = mcon({tr{1},tr{2},tr{3},tr{4}},...
%     {[1,5,2,5],[2,6,3,6],[3,7,4,7],[4,8,1,8]});
% Z1 = Z1*prod(tf);
% disp((Z0-Z1)/Z0);

% z1 = fn_original(lambda,B,M,N);
z1 = trace((m1*m2)^(M/2));
disp(z1);

D = 40;
tr = tensornet(1,:);
m_factor = zeros(M-1,N);
tot_factor = 1;
for k = 2:M
   [tr,tf] = fn_add_line2(tr,tensornet(k,:),D);
   m_factor(k-1,:) = tf;
   tot_factor = tot_factor*prod(tf);
end


tv1 = 1:N;
tv2 = (N+1):(2*N);
tv3 = [2:N,1];
tv4 = (N+1):(2*N);
tm = [tv1;tv2;tv3;tv4];
tm = transpose(tm);
cind = cell(1,N);
for k=1:N
   cind{k} = tm(k,:); 
end
z2 =  mcon(tr,cind);
z2 = z2*tot_factor;
disp(['z2=',num2str(z2)]);
disp(['z1=',num2str(z1)]);
disp(['computatioal error=',num2str((z2-z1)/z1)]);
disp(['max svd error=',num2str(max_error)]);

% %%%%%%%%% test code %%%%%%%%%%%%%%%%%%%%%%%%%%
% m_configuration = sign(rand(M,N)-0.5);
% m_op = m_configuration;
% x = floor(rand()*N) + 1;
% y = floor(rand()*M) + 1;
% [res,num_circle,m_pearl] = fn_factor(x,y,m_op,B);
% disp(m_op);
% fprintf('x = %d ; y = %d \n',x,y);
% fprintf('res = %f\n',res);
% fprintf('number of circles = %d\n',num_circle);
% fprintf('pearl matrix: \n');
% disp(m_pearl);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% z1 = fn_original(lambda,B,M,N);
% z2 = fn_original(lambda+1e-6,B,M,N);
% partial_theta_exact = (1-lambda)^2*(z2-z1)/z1/1e-6;
% 
% z1 = fn_original(lambda,B,M,N);
% z2 = fn_original(lambda,B+1e-6,M,N);
% partial_B_exact = (z2-z1)/z1/1e-6;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cTr,v_factor] = fn_add_line(cT,cA,D)

N = length(cT);
cTr = cell(1,N);
v_factor = zeros(1,N);
T = cT{1}; A = cA{1};
TA = ncon({T,A},{[-1,-2,-3,1],[-6,1,-4,-5]});
TA = permute(TA,[1,6,2,3,4,5]);
ts = size(TA);
TA = reshape(TA,[ts(1)*ts(2),ts(3),ts(4),ts(5),ts(6)]);

T_right = TA;

for n=2:N
   T_left = T_right; 
   T = cT{n}; A = cA{n};
   TA = ncon({T,A},{[-1,-2,-3,1],[-6,1,-4,-5]});
   T_right = TA;
   [T_left,T_right,factor] = fn_update_LR(T_left,T_right,[3,4],[1,6],D);
   T_left = permute(T_left,[1,2,4,3]);
   cTr{n-1} = T_left;
   v_factor(n-1) = factor;
end
T_left = T_right;
ts = size(T_left);
T_left = reshape(T_left,[ts(1),ts(2),ts(3)*ts(4),ts(5)]);
T_right = cTr{1};
[T_left,T_right,factor] = fn_update_LR(T_left,T_right,[3],[1],D);
T_left = permute(T_left,[1,2,4,3]);
cTr{N} = T_left;
cTr{1} = T_right;
v_factor(N) = factor;
end

function [TL_new,TR_new,factor] = fn_update_LR(TL,TR,vL,vR,D)

tsL = size(TL);
tsR = size(TR);
indL = 1:length(tsL);
indR = 1:length(tsR);
vL_L = indL; vL_L(vL) = [];
vR_R = indR; vR_R(vR) = [];
TL = permute(TL,[vL_L,vL]);
TR = permute(TR,[vR,vR_R]);
mL = reshape(TL,prod(tsL(vL_L)),prod(tsL(vL)));
mR = reshape(TR,prod(tsR(vR)),prod(tsR(vR_R)));
truncate_D = min(D,prod(tsR(vR)));
[uLm,sLm,vLm] = svd(mL,'econ'); vLm = vLm';
[uRm,sRm,vRm] = svd(mR,'econ'); vRm = vRm';
uLm = uLm(:,1:truncate_D);
%
tv = diag(sLm);
tv2 = zeros(1,length(tv));
tv2(1:truncate_D) = tv(1:truncate_D);
disp(sqrt(norm(tv)^2-norm(tv2)^2)/norm(tv));
%
sLm = sLm(1:truncate_D,1:truncate_D);
vLm = vLm(1:truncate_D,:);
uRm = uRm(:,1:truncate_D);
%
tv = diag(sRm);
tv2 = zeros(1,length(tv));
tv2(1:truncate_D) = tv(1:truncate_D);
disp(sqrt(norm(tv)^2-norm(tv2)^2)/norm(tv));
%
sRm = sRm(1:truncate_D,1:truncate_D);
vRm = vRm(1:truncate_D,:);
factor = norm(diag(sLm))*norm(diag(sRm));
mL = uLm;
mR = sLm*vLm*uRm*sRm*vRm;
mR = mR/factor;
TL_new = reshape(mL,[tsL(vL_L),truncate_D]);
TR_new = reshape(mR,[truncate_D,tsR(vR_R)]);

end


function [cTr,v_factor] = fn_add_line2(cT,cA,D)

N = length(cT);
cTr = cell(1,N);
v_factor = zeros(1,N);
T = cT{1}; A = cA{1};
TA = ncon({T,A},{[-1,-2,-3,1],[-6,1,-4,-5]});
TA = permute(TA,[1,6,2,3,4,5]);
ts = size(TA);
TA = reshape(TA,[ts(1)*ts(2),ts(3),ts(4),ts(5),ts(6)]);

T_right = TA;

for n=2:N
   T_left = T_right; 
   T = cT{n}; A = cA{n};
   TA = ncon({T,A},{[-1,-2,-3,1],[-6,1,-4,-5]});
   T_right = TA;
   [T_left,T_right,factor] ...
       = fn_update_LR2(T_left,T_right,[3,4],[1,6],D);
   T_left = permute(T_left,[1,2,4,3]);
   cTr{n-1} = T_left;
   v_factor(n-1) = factor;
end
T_left = T_right;
ts = size(T_left);
T_left = reshape(T_left,[ts(1),ts(2),ts(3)*ts(4),ts(5)]);
T_right = cTr{1};
[T_left,T_right,factor] = fn_update_LR2(T_left,T_right,[3],[1],D);
T_left = permute(T_left,[1,2,4,3]);
cTr{N} = T_left;
cTr{1} = T_right;
v_factor(N) = factor;
end


function [TL_new,TR_new,factor] = fn_update_LR2(TL,TR,vL,vR,D)
global max_error;
tsL = size(TL);
tsR = size(TR);
indL = 1:length(tsL);
indR = 1:length(tsR);
vL_L = indL; vL_L(vL) = [];
vR_R = indR; vR_R(vR) = [];
TL = permute(TL,[vL_L,vL]);
TR = permute(TR,[vR,vR_R]);
mL = reshape(TL,prod(tsL(vL_L)),prod(tsL(vL)));
mR = reshape(TR,prod(tsR(vR)),prod(tsR(vR_R)));
matrix = mL*mR;
truncate_D = min([D,prod(tsR(vR_R)),prod(tsL(vL_L))]);

[u,s,v] = svd(matrix,'econ'); v = v';
ut = u(:,1:truncate_D);
%
tv = diag(s);
tv2 = zeros(1,length(tv));
tv2(1:truncate_D) = tv(1:truncate_D);
error_here = sqrt(norm(tv)^2-norm(tv2)^2)/norm(tv);
disp(['truncate error=',num2str(error_here)]);
max_error = max(max_error,error_here);
%
st = s(1:truncate_D,1:truncate_D);
vt = v(1:truncate_D,:);
factor = norm(diag(st));
mL = ut;
mR = st*vt;
mR = mR/factor;
TL_new = reshape(mL,[tsL(vL_L),truncate_D]);
TR_new = reshape(mR,[truncate_D,tsR(vR_R)]);

end


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