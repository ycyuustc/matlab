close(figure(1));
figure(1);
hold on;

Num_k = 400;

M = 40;
Num = 2*M;
gamma = 3;

vt1 = linspace(-pi,0,Num_k);

for i=1:Num_k

t1 = vt1(i);

a = t1 - gamma/2;
b = t1 + gamma/2;
c = 1;
V = 1;

% a = 1.5;
% b = 0.5;
% c = 1;

tv1 = reshape([a*ones(1,M);c*ones(1,M)],[1,Num]);
tv2 = reshape([b*ones(1,M);c*ones(1,M)],[1,Num]);


m = diag(tv1(1:end),-1) + diag(tv2(1:end),1);
% m(1,Num) = c; m(Num,1) = c;  m(1,1) = V;
% disp('H='); disp(num2str(m));

v_eig = eigs(m,Num+1);

% tv =  v(:,Num)/v(2,8);
% lambda = d(Num,Num);
% 
% A = b;
% B = a*b + 1 - lambda^2;
% C = a;
% 
% disp(A*tv(8)+B*tv(6)+C*tv(4));
% 
% Del = B^2-4*A*C;
% if Del>=0
%     x1 =1/2/A*(-B+sqrt(Del));
%     x2 =1/2/A*(-B-sqrt(Del));
% else
%     x1 =1/2/A*(-B+1j*sqrt(abs(Del)));
%     x2 =1/2/A*(-B-1j*sqrt(abs(Del)));
% end
figure(1);
hold on;
plot(t1*ones(1,Num+1),abs(v_eig),'b*');


end