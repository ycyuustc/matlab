m2id=eye(2);
m2z=[1,0;0,-1];

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

lambda = 0.1;

p1 = 1/(1-lambda);
p2 = lambda/(lambda-1);

T1 = p1*m4p + p2*m4id;
T2 = p1*m4p2 + p2*m4id;

t = reshape(T1,4,4);
[u,s,v] = svd(t);
disp(diag(s));
u1 = u*sqrt(s);
v1 = sqrt(s)*v';
A1 = reshape(u1,[2,2,4]);
D1 = reshape(v1,[4,2,2]);
D1 = permute(D1,[2,3,1]);

t = permute(T1,[1,4,2,3]);
t = reshape(t,4,4);

[u,s,v] = svd(t);
disp(diag(s));
u1 = u*sqrt(s);
v1 = sqrt(s)*v';
B1 = reshape(v1,[4,2,2]);
B1 = permute(B1,[2,3,1]);
C1 = reshape(u1,[2,2,4]);
C1 = permute(C1,[2,1,3]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = reshape(T2,4,4);
[u,s,v] = svd(t);
disp(diag(s));
u1 = u*sqrt(s);
v1 = sqrt(s)*v';
A2 = reshape(u1,[2,2,4]);
D2 = reshape(v1,[4,2,2]);
D2 = permute(D2,[2,3,1]);

t = permute(T2,[1,4,2,3]);
t = reshape(t,4,4);
[u,s,v] = svd(t);
disp(diag(s));
u1 = u*sqrt(s);
v1 = sqrt(s)*v';
B2 = reshape(v1,[4,2,2]);
B2 = permute(B2,[2,3,1]);
C2 = reshape(u1,[2,2,4]);
C2 = permute(C2,[2,1,3]);

%%% check:
fn_c=@(x1,x2,x3,x4,x5,x6) fn_contract(x1,x2,x3,x4,x5,x6);
t = fn_c(A1,3,3,D1,3,3);
disp(norm(reshape((t - T1),4,4)));
t = fn_c(B1,3,3,C1,3,3);
t = permute(t,[4,1,2,3]);
disp(norm(reshape(t-T1,4,4)));

t = fn_c(A2,3,3,D2,3,3);
disp(norm(reshape((t - T2),4,4)));
t = fn_c(B2,3,3,C2,3,3);
t = permute(t,[4,1,2,3]);
disp(norm(reshape(t-T2,4,4)));

%%%%%%%%  end check  %%%%%%%%%%%%%%%%%55

% calculate the big T
A = A1; B = B1; C = C2; D = D2;
T = D;
T = fn_c(T,3,1,C,3,2);
T = fn_c(T,4,3,A,3,2);
T = fn_c(T,5,[1,4],B,3,[1,2]);
bT1 = T;

A = A2; B = B2; C = C1; D = D1;
T = D;
T = fn_c(T,3,1,C,3,2);
T = fn_c(T,4,3,A,3,2);
T = fn_c(T,5,[1,4],B,3,[1,2]);
bT2 = T;

t = reshape(bT1,16,16);
[u,s,v] = svd(t);
u = u(:,1:4);
s = s(1:4,1:4);
v = v(:,1:4);
u1 = u*sqrt(s);
v1 = sqrt(s)*v';

bA1 = reshape(u1,[4,4,4]);
bD1 = reshape(v1,[4,4,4]);
bD1 = permute(bD1,[2,3,1]);

t = permute(bT2,[1,4,2,3]);
t = reshape(t,16,16);
[u,s,v] = svd(t);
u = u(:,1:4);
s = s(1:4,1:4);
v = v(:,1:4);
u1 = u*sqrt(s);
v1 = sqrt(s)*v';
bB2 = reshape(v1,[4,4,4]);
bB2 = permute(bB2,[2,3,1]);
bC2 = reshape(u1,[4,4,4]);
bC2 = permute(bC2,[2,1,3]);

%%%%%   check %%%%%%%%%%%
t = fn_c(bB2,3,3,bC2,3,3);
t = permute(t,[4,1,2,3]);
disp(norm(reshape(t-bT2,16,16)));
t = fn_c(bA1,3,3,bD1,3,3);
disp(norm(reshape(t-bT1,16,16)));
%%%%%  end of the check %%%%%%%%%%%%%

%%%%%  calculate the big big T  %%%%%
A = bA1; B = bB2; C = bC2; D = bD1;
T = D;
T = fn_c(T,3,1,C,3,2);
T = fn_c(T,4,3,A,3,2);
T = fn_c(T,5,[1,4],B,3,[1,2]);
BT = T;

truncate = 8;





