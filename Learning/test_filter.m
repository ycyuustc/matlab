m4p=zeros(2,2,2,2);
m4p(1,1,1,1)=1;
m4p(1,2,2,1)=1;
m4p(2,1,1,2)=1;
m4p(2,2,2,2)=1;

m4p2 = permute(m4p,[1,4,3,2]);
theta = 0.01;
R1 = (1-theta)*m4p + theta*m4p2;
R2 = (1-theta)*m4p2 + theta*m4p;

tm1 = permute(R1,[1,4,3,2]);
tm1 = reshape(tm1,4,4);
[u,s,v] = svd(tm1,'econ'); v = v';
u = u*sqrt(s);
v = sqrt(s)*v;
C = reshape(u,[2,2,4]);
A = reshape(v,[4,2,2]);
A = permute(A,[2,3,1]);

tm2 = reshape(R2,4,4);
[u,s,v] = svd(tm2,'econ'); v = v';
u = u*sqrt(s);
v = sqrt(s)*v;
B = reshape(u,[2,2,4]);
B = permute(B,[2,1,3]);
D = reshape(v,[4,2,2]);
D = permute(D,[3,2,1]);

% t1 = ncon({A,C},{[-3,-2,1],[-1,-4,1]});
% t2 = ncon({B,D},{[-2,-1,1],[-4,-3,1]});
% T1 = ncon({v,u},{[-1,1,-3],[1,-2,-4]});
% T2 = ncon({v,u},{[});

T1 = ncon({A,C},{[1,-1,-3],[1,-2,-4]});
T2 = ncon({A,C},{[-1,1,-4],[-2,1,-3]});

ts = size(T1);
tm = reshape(T1,ts(1)*ts(2),ts(3)*ts(4));
[u,s,v] = svd(tm,'econ'); v = v';
u = u*sqrt(s);
v = sqrt(s)*v;
t1 = reshape(u,[ts(1),ts(2),numel(u)/ts(1)/ts(2)]);

ts = size(T2);
tm = reshape(T2,ts(1)*ts(2),ts(3)*ts(4));
[u,s,v] = svd(tm,'econ'); v = v';
u = u*sqrt(s);
v = sqrt(s)*v;
t2 = reshape(u,[ts(1),ts(2),numel(u)/ts(1)/ts(2)]);

t12 = ncon({t1,t2},{[-1,1,-3],[1,-2,-4]});


T1T1 = ncon({T1,T1},{[-1,-3,-5,1],[-2,-4,1,-6]});
ts2 = size(T1T1); ts = [ts2,ones(6-length(ts2))];
T1T1 = reshape(T1T1,[ts(1)*ts(2),ts(3)*ts(4),ts(5),ts(6)]);
tm = reshape(T1T1,[ts(1)*ts(2)*ts(3)*ts(4),ts(5)*ts(6)]);
[u,s,v] = svd(tm,'econ'); v = v';

T1T2 = ncon({T1,T2},{[-1,1,-3,-5],[1,-2,-4,-6]});
ts2 = size(T1T2); ts = [ts2,ones(6-length(ts2))];
T1T2 = reshape(T1T2,[ts(1),ts(2),ts(3)*ts(4),ts(5)*ts(6)]);
tm = reshape(T1T2,[ts(1)*ts(2),ts(3)*ts(4)*ts(5)*ts(6)]);
[u,s,v] = svd(tm,'econ'); v = v';
u = u*sqrt(s);
v = sqrt(s)*v;
t = reshape(v,[numel(v)/ts(3)/ts(4)/ts(5)/ts(6),ts(3)*ts(4),ts(5)*ts(6)]);
t = permute(t,[2,3,1]);
[t,r] = fn_qr(eye(16),t);


A = rand(2,2,2)+rand(2,2,2)*1j;
B = rand(2,2,2)+rand(2,2,2)*1j;
lambda = diag(rand(1,2));
theta = diag(rand(1,2));

[A,lambda,a] = fn_TEBD_canonical(A,lambda);
[B,theta,b] = fn_TEBD_canonical(B,theta);
disp(ncon({A,conj(A),lambda,conj(lambda)},{[-1,2,1],[-2,3,1],[2,4],[3,4]}));
disp(ncon({B,conj(B),theta,conj(theta)},{[-1,2,1],[-2,3,1],[2,4],[3,4]}));
disp(ncon({A,conj(A),lambda,conj(lambda)},{[2,-1,4],[3,-2,4],[1,2],[1,3]}));
disp(ncon({B,conj(B),theta,conj(theta)},{[2,-1,4],[3,-2,4],[1,2],[1,3]}));

[A1,lambda1,B1,theta1,sR,sL,etaR,etaL] = fn_TEBD_canonical_bi(A,lambda,B,theta);
disp('haha');
disp(ncon({A1,conj(A1),lambda1,conj(lambda1)},{[-1,2,1],[-2,3,1],[2,4],[3,4]}));
disp(ncon({B1,conj(B1),theta1,conj(theta1)},{[-1,2,1],[-2,3,1],[2,4],[3,4]}));
disp(ncon({A1,conj(A1),lambda1,conj(lambda1)},{[2,-1,4],[3,-2,4],[1,2],[1,3]}));
disp(ncon({B1,conj(B1),theta1,conj(theta1)},{[2,-1,4],[3,-2,4],[1,2],[1,3]}));

disp(ncon({A1,conj(B1),lambda1,conj(theta1),sR},...
    {[-1,2,1],[-2,3,1],[2,4],[3,5],[4,5]}));
disp(diag(ncon({A1,conj(B1),lambda1,conj(theta1),sR},...
    {[-1,2,1],[-2,3,1],[2,4],[3,5],[4,5]}))./diag(sR));
disp(etaR);
disp(ncon({A1,conj(B1),lambda1,conj(theta1),sL},...
    {[2,-1,5],[4,-2,5],[1,2],[3,4],[3,1]}));
disp(diag(ncon({A1,conj(B1),lambda1,conj(theta1),sL},...
    {[2,-1,5],[4,-2,5],[1,2],[3,4],[3,1]}))./diag(sL));
disp(etaL);


hA = ncon({lambda,conj(lambda),A,conj(A),lambda,conj(lambda)},...
    {[1,2],[1,3],[2,4,-2],[3,5,-1],[4,6],[5,6]});
disp(hA);
hA1 = ncon({lambda1,conj(lambda1),A1,conj(A1),lambda1,conj(lambda1)},...
    {[1,2],[1,3],[2,4,-2],[3,5,-1],[4,6],[5,6]});
disp(hA1);
hB = ncon({theta,conj(theta),B,conj(B),theta,conj(theta)},...
    {[1,2],[1,3],[2,4,-2],[3,5,-1],[4,6],[5,6]});
disp(hB);
hB1 = ncon({theta1,conj(theta1),B1,conj(B1),theta1,conj(theta1)},...
    {[1,2],[1,3],[2,4,-2],[3,5,-1],[4,6],[5,6]});
disp(hB1);

% while 1
%     
%    [t,r] = fn_qr(r,t); 
%    disp(r);
%     
% end

function [Q,R] = fn_qr(T1,T2)

t = mcon({T1,T2},{[-1,1],[1,-2,-3]});
vs = size(t);
t = reshape(permute(t,[1,3,2]),vs(1)*vs(3),vs(2));
[Q,R]=qr(t);
Q = Q(:,1:vs(2));
R = R(1:vs(2),:);

v_diag = sign(diag(R));
tran_m = diag(v_diag);
Q = Q*tran_m;
R = tran_m*R;

Q = reshape(Q,[vs(1),vs(3),vs(2)]);
Q = permute(Q,[1,3,2]);
R = reshape(R,vs(2),vs(2));

end
