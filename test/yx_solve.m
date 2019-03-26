Num = 16;
V = 1;
g = 0.01;

a = exp(g*Num);

matrix = diag(-ones(1,2*Num+1),-1);
va = zeros(1,2*Num+2);
va(1) = -V;
va(2) = -1;
va(Num) = -(a^2+1)/a;
va(Num+2) = (a^2+1)/a;
va(2*Num) = 1;
va(2*Num+1) = V;
va(2*Num+2) = -1;
va = va.*power(-1,1:(2*Num+2));
matrix(1,:) = va;

ve = eig(matrix);
fprintf('The solutions are:\n');
disp(ve);

%%%% test the solution
for k = 1:length(ve)
    p = ve(k);
    q = 1/p;
    LHS = 1/p + V/(a*p^Num-1);
    RHS = 1/q + V/(a*q^Num-1);
    flag = abs(LHS-RHS);
    fprintf('The checking error for the %dth root is %f\n',k,flag);
end