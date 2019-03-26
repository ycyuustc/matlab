close(figure(1));
figure(1);
% hold on;

x = 1.1*pi*rand();
K = 2*rand();
delta_x = 0.01;
disp(fn_fj2(3,K,x));
disp((fn_fj2(3,K,x+0.00001)-fn_fj2(3,K,x-0.00001))/0.00002);
disp(fn_fj2_dx(3,K,x));

disp((fn_fj2(3,K,x+delta_x)+fn_fj2(3,K,x-delta_x)-2*fn_fj2(3,K,x))/delta_x^2);
disp(fn_fj2_dxdx(3,K,x));

disp(fn_fj(3,1,x));
disp(fn_fj2(3,1,x));

M = 20;
Num = M*2;

gamma = 3;

v_t1 = linspace(sqrt(gamma^2/4+1),gamma/2+0.001,100);
vx0 = exp(1j*2*pi/(2*M+1)*(1:1:(2*M)));
figure(1);
plot(real(vx0),imag(vx0),'r*');
for i=2:100
    
    t1 = v_t1(i);
    t2 = 1;
    a = t1 - gamma/2;
    b = t1 + gamma/2;
    K = a*b;
    K = 1/sqrt(K);
    
    vx0 = fn_solve22(M,K,vx0);
    disp(vx0);
    figure(1);
    plot(real(vx0),imag(vx0),'r*');
    disp(i);
    
end

res_check = zeros(1,2*M);
for i = 1:(2*M)
    x = vx0(i);
    res_check(i) = fn_fj(M,K,x);
end
disp(res_check);

disp(sqrt(K+1+res+K./res));

v_energy = sqrt(K+1+res+K./res);
disp(v_energy);



% a =
b = 1/a;
c = 1;

tv1 = reshape([a*ones(1,M);c*ones(1,M)],[1,Num]);
tv2 = reshape([b*ones(1,M);c*ones(1,M)],[1,Num]);


m = diag(tv1(1:end-1),-1) + diag(tv2(1:end-1),1);
% m(1,Num) = c; m(Num,1) = c;  m(1,1) = V;
disp('H='); disp(num2str(m));

[v,d] = eig(m);

tv =  v(:,Num)/v(2,Num);
lambda = d(Num,Num);

A = b;
B = a*b + 1 - lambda^2;
C = a;

% disp(A*tv(8)+B*tv(6)+C*tv(4));

Del = B^2-4*A*C;
if Del>=0
    x1 =1/2/A*(-B+sqrt(Del));
    x2 =1/2/A*(-B-sqrt(Del));
else
    x1 =1/2/A*(-B+1j*sqrt(abs(Del)));
    x2 =1/2/A*(-B-1j*sqrt(abs(Del)));
end

K = 0.2;
M=8;
res = fn_solve(M,K);
res_check = zeros(1,M);
for i = 1:M
    x = res(i);
    res_check(i) = fn_fj(M,K,x);
end
disp(res_check);
disp(sqrt(K+1+res+K./res));

v_energy = sqrt(K+1+res+K./res);
disp(v_energy);

function vx_res = fn_solve22(N,a,vx0)

% vx0 = exp(1j*2*pi/(2*N+1)*(1:1:(2*N)));
vx_res = zeros(1,2*N);
for n = 1:(2*N)   
    x0 = vx0(n); 
    if real(x0)<-1.2 && abs(imag(x0))<1e-8
        vx_res(n) = x0;
    else
    x0_pre = x0;   
    
    x0 = fn_step2(N,a,x0);
    
    if abs(imag(x0_pre))>1e-6 && abs(imag(x0))<1e-6
        
        if imag(x0_pre)>0
            x0 = 0;
            x0 = fn_step2(N,a,x0);
        else
            x0 = x0 + 0.2;
            x0 = fn_step2(N,a,x0);
        end
        
    end
    
    vx_res(n) = x0;
    end
    
end

end

function vx_res = fn_solve2(N,K,vx0)

% vx0 = exp(1j*2*pi/(2*N+1)*(1:1:(2*N)));
vx_res = zeros(1,2*N);
for n = 1:(2*N)   
    x0 = vx0(n);  
    x0_pre = x0;   
    
    x0 = fn_step(N,K,x0);
    
    if abs(imag(x0_pre))>1e-6 && abs(imag(x0))<1e-6
        
        if imag(x0_pre)>0
            x0 = 0;
            x0 = fn_step(N,K,x0);
        else
            x0 = x0 + 0.2;
            x0 = fn_step(N,K,x0);
        end
        
    end
    
    vx_res(n) = x0;
    
end

end

function x0 = fn_step2(N,a,x0)

A = 1/2*fn_fj2_dxdx(N,a,x0);
B = fn_fj2_dx(N,a,x0);
C = fn_fj2(N,a,x0);
while abs(C)>1e-8
    if abs(C/B)>1e-8
        delta_x = -C/B;
    else
        x1 = 1/2/A*(-B+sqrt(B^2-4*A*C));
        x2 = 1/2/A*(-B-sqrt(B^2-4*A*C));
        C1 = fn_fj2(N,a,x0+x1);
        C2 = fn_fj2(N,a,x0+x2);
        if abs(abs(C1)<abs(C2))
            delta_x = x1;
        else
            delta_x = x2;
        end
    end
    
    x0 = x0 + delta_x;
    A = 1/2*fn_fj2_dxdx(N,a,x0);
    B = fn_fj2_dx(N,a,x0);
    C = fn_fj2(N,a,x0);
    disp(abs(C));
    
end

end


function x0 = fn_step(N,K,x0)

A = 1/2*fn_fj_dxdx(N,K,x0);
B = fn_fj_dx(N,K,x0);
C = fn_fj(N,K,x0);
while abs(C)>1e-8
    if abs(B)>1e-8
        delta_x = -C/B;
    else
        x1 = 1/2/A*(-B+sqrt(B^2-4*A*C));
        x2 = 1/2/A*(-B-sqrt(B^2-4*A*C));
        if abs(x1<x2)
            delta_x = x1;
        else
            delta_x = x2;
        end
    end
    
    x0 = x0 + delta_x;
    A = 1/2*fn_fj_dxdx(N,K,x0);
    B = fn_fj_dx(N,K,x0);
    C = fn_fj(N,K,x0);
    
end

end


function vx_res = fn_solve(N,K)
Num = 100;
vK = linspace(1,K,Num);

vx0 = exp(1j*2*pi/(2*N+1)*(1:1:(2*N)));
vx_res = zeros(1,2*N);
for n = 1:(2*N)
    
    x0 = vx0(n);
    
    for ik = 2:Num
        K = vK(ik);
        A = 1/2*fn_fj_dxdx(N,K,x0);
        B = fn_fj_dx(N,K,x0);
        C = fn_fj(N,K,x0);
        while abs(C)>1e-8
            if abs(B)>1e-5
                delta_x = -C/B;
            else
                x1 = 1/2/A*(-B+sqrt(B^2-4*A*C));
                x2 = 1/2/A*(-B-sqrt(B^2-4*A*C));
                if abs(x1<x2)
                    delta_x = x1;
                else
                    delta_x = x2;
                end
            end
            
            x0 = x0 + delta_x;
            A = 1/2*fn_fj_dxdx(N,K,x0);
            B = fn_fj_dx(N,K,x0);
            C = fn_fj(N,K,x0);
            
        end
        
    end
    
    vx_res(n) = x0;
    
end

end

function res = fn_fj2(N,a,x)

vx = power(x,(2*N:-1:0));
va = reshape([ones(1,N+1);a*ones(1,N+1)],[1,(N+1)*2]);
va = va(1:end-1);
res = sum(vx.*va);

end

function res = fn_fj(N,K,x)

vx = power(x,(2*N:-1:0));
vK = reshape([power(K,0:1:N);power(K,0:1:N)],[1,(N+1)*2]);
vK = vK(1:end-1);
res = sum(vx.*vK);

end

function res = fn_fj2_dx(N,a,x)

vx = power(x,((2*N-1):-1:0));
va = reshape([ones(1,N);a*ones(1,N)],[1,(N)*2]);
vN = (2*N):-1:1;

res = sum(vx.*va.*vN);

end


function res = fn_fj_dx(N,K,x)

vx = power(x,((2*N-1):-1:0));
vK = reshape([power(K,0:1:(N-1));power(K,0:1:(N-1))],[1,N*2]);
vN = (2*N):-1:1;

res = sum(vx.*vK.*vN);

end


function res = fn_fj2_dxdx(N,a,x)

vx = power(x,((2*N-2):-1:0));
va = reshape([ones(1,N);a*ones(1,N)],[1,(N)*2]);
va = va(1:end-1);
vN1 = (2*N):-1:2;
vN2 = (2*N-1):-1:1;

res = sum(vx.*va.*vN1.*vN2);

end


function res = fn_fj_dxdx(N,K,x)

vx = power(x,((2*N-2):-1:0));
vK = reshape([power(K,0:1:(N-1));power(K,0:1:(N-1))],[1,N*2]);
vK = vK(1:end-1);
vN1 = (2*N):-1:2;
vN2 = (2*N-1):-1:1;

res = sum(vx.*vK.*vN1.*vN2);

end


