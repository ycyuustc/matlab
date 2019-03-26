
L=1;
% k=zeros(9,8);

% load('vgg4.mat')
% k=vg;

% load('vg3.mat')
% k=vg;

% load('ve1_9_c_k.mat')
% k=ve1;
Num_p = 100;
v_p = linspace(-8*pi,8*pi,Num_p);
v_res = zeros(1,Num_p);
vk1 = kk_t12(:,1);
vk2 = kk_t12(:,2);

% for i=1:Num_p
%    momentum = v_p(i);
%    v_res(i) = fn_int(momentum);
%    disp(i);
%    disp(v_res(i));
% end
% 
% figure(1);
% plot(v_p,real(v_res),'bo');
% figure(2);
% plot(v_p,imag(v_res),'ro');

L = 1;
k=[1.0*pi 5.0/3*pi];
k10=k(1,1);k20=k(1,2);

k11=(k10+3*k20)/2; k21=(k10-k20)/2; k12=(k10-3*k20)/2; k22=(k10+k20)/2;
%A??????????A0+??A1+??A2+??A0-??A1-??A2-
A=[1,1,0,0,0,0;0,0,0,-exp(1i*(k20-k10)*L),-exp(1i*(k21-k11)*L),0;-exp(-1i*k20*L),...
    0,1,0,0,0;0,0,0,exp(-1i*k10*L),0,-exp(1i*(k22-k12)*L);0,-exp(-1i*k21*L),-exp(1i*k12*L),...
    0,0,0;0,0,0,0,exp(-1i*k11*L),exp(1i*k22*L)];
B1=[A;1,0,0,-1,0,0];
bs=null(B1);
%A1j1=Aj+(+,+) A1j2=Aj+(+,-) A1j3=Aj+(-,+) A1j4=Aj+(-,-) A2j1=Aj-(+,+) A2j2=Aj-(+,-) A2j3=Aj-(-,+) A2j4=Aj-(-,-)
A101=bs(1);A111=bs(2);A121=bs(3);A201=bs(4);A211=bs(5);A221=bs(6);
A102=-exp(-1i*k20*L)*bs(1);A112=-exp(-1i*k21*L)*bs(2);A122=-exp(-1i*k22*L)*bs(3);
A202=-exp(1i*k20*L)*bs(4);A212=-exp(1i*k21*L)*bs(5);A222=-exp(1i*k22*L)*bs(6);
A103=-exp(1i*k10*L)*bs(1);A113=-exp(1i*k11*L)*bs(2);A123=-exp(1i*k12*L)*bs(3);
A203=-exp(-1i*k10*L)*bs(4);A213=-exp(-1i*k11*L)*bs(5);A223=-exp(-1i*k12*L)*bs(6);
A104=exp(1i*k10*L-1i*k20*L)*bs(1);A114=exp(1i*k11*L-1i*k21*L)*bs(2);
A124=exp(1i*k12*L-1i*k22*L)*bs(3);A204=exp(1i*k20*L-1i*k10*L)*bs(4);
A214=exp(1i*k21*L-1i*k11*L)*bs(5);A224=exp(1i*k22*L-1i*k12*L)*bs(6);

Num = 501;
vx1 = linspace(-0.5,0.5,Num);
vx2 = linspace(-0.5,0.5,Num);
[x1,x2] = meshgrid(vx1,vx2);
f=(A101*exp(1i*k10*x1+1i*k20*x2)+A102*exp(1i*k10*x1-1i*k20*x2)+A103*exp(-1i*k10*x1+1i*k20*x2)+A104*exp(-1i*k10*x1-1i*k20*x2)+...
    A111*exp(1i*k11*x1+1i*k21*x2)+A112*exp(1i*k11*x1-1i*k21*x2)+A113*exp(-1i*k11*x1+1i*k21*x2)+A114*exp(-1i*k11*x1-1i*k21*x2)+...
    A121*exp(1i*k12*x1+1i*k22*x2)+A122*exp(1i*k12*x1-1i*k22*x2)+A123*exp(-1i*k12*x1+1i*k22*x2)+A124*exp(-1i*k12*x1-1i*k22*x2)).*(x2-x1<0)+...
    (A201*exp(1i*k10*x1+1i*k20*x2)+A202*exp(1i*k10*x1-1i*k20*x2)+A203*exp(-1i*k10*x1+1i*k20*x2)+A204*exp(-1i*k10*x1-1i*k20*x2)+...
    A211*exp(1i*k11*x1+1i*k21*x2)+A212*exp(1i*k11*x1-1i*k21*x2)+A213*exp(-1i*k11*x1+1i*k21*x2)+A214*exp(-1i*k11*x1-1i*k21*x2)+...
    A221*exp(1i*k12*x1+1i*k22*x2)+A222*exp(1i*k12*x1-1i*k22*x2)+A223*exp(-1i*k12*x1+1i*k22*x2)+A224*exp(-1i*k12*x1-1i*k22*x2)).*(x2-x1>=0);

delta_x = 1/(Num-1);
m_nt = ones(Num,Num);
m_nt(1,1) = 0.25;
m_nt(1,Num) = 0.25;
m_nt(Num,1) = 0.25;
m_nt(Num,Num) = 0.25;
m_nt(1,2:(Num-1)) = 0.5;
m_nt(Num,2:(Num-1)) = 0.5;
m_nt(2:(Num-1),1) = 0.5;
m_nt(2:(Num-1),Num) = 0.5;

m_nt = m_nt * delta_x^2;
figure(2);
gca = pcolor(x1,x2,imag(f));
set(gca,'linestyle','none');
% 
% disp(interp2(x1,x2,conj(f).*f,0,0.3,'Nearest'));
load('E:\matlab\test\A.mat');
% disp(fn_density(0,0.3,A,vk1,vk2));

k1 = rand()-1/2;
k2 = rand()-1/2;
int_1 = sum(sum(m_nt.*conj(f).*f.*exp(1j*k1*x1-1j*k2*x2)));
disp(int_1);
disp(fn_density_k(k1,k2,A,vk1,vk2));

k1 = 2*pi*8;
k2 = 2*pi*3;

int_1 = sum(sum(m_nt.*f.*exp(1j*k1*x1-1j*k2*x2)));
disp(abs(int_1));

k1 = -k1;
k2 = -k2;

int_2 = sum(sum(m_nt.*f.*exp(1j*k1*x1-1j*k2*x2)));
disp(abs(int_2));


function res = fn_density_k(p1,p2,A,vk1,vk2)

v_sigma = [1,-1];
my_sum = 0;

m_permutation = perms([2,1]);

for im = 1:2
    
v_sort = m_permutation(im,:);
    
for i1 = 1:3
    for i2 = 1:3
        for si11 = 1:2
            for si21 = 1:2
                for si12 = 1:2
                    for si22 = 1:2
                        
                        vx = zeros(1,2);
                        vx(v_sort) = [1,2];
                        
                        if vx(1)>vx(2)
                            r1 = 1; r2 = 1;
                        else
                            r1 = 2; r2 = 2;
                        end
                                           
                        sigma11 = v_sigma(si11);
                        sigma12 = v_sigma(si12);
                        sigma21 = v_sigma(si21);
                        sigma22 = v_sigma(si22);
                        
                        k1 = -sigma11*vk1(i1) + sigma12*vk1(i2) + p1;
                        k2 = -sigma21*vk2(i1) + sigma22*vk2(i2) - p2;
                        vk = [k1,k2];
                        vk = vk(v_sort);
                        k1 = vk(1); k2 = vk(2);
                        
                        my_sum = my_sum + ...
                            A(r2,i2,si12,si22)*...
                            conj(A(r1,i1,si11,si21))*...
                            fn_g([k1+k2,k2,0])*...
                            exp(-1j/2*(k1+k2));
                        
                    end
                end
            end
        end
    end
end

end

res = my_sum;


end
































function res = fn_density(x1,x2,A,vk1,vk2)
v_sigma = [1,-1];
my_sum = 0;
for i1 = 1:3
    for i2 = 1:3
        for si11 = 1:2
            for si21 = 1:2
                for si12 = 1:2
                    for si22 = 1:2
                        
                        if x1>x2
                            r1 = 1;
                        else
                            r1 = 2;
                        end
                        
                        if x1>x2
                            r2 = 1;
                        else
                            r2 = 2;
                        end
                        
                        sigma11 = v_sigma(si11);
                        sigma12 = v_sigma(si12);
                        sigma21 = v_sigma(si21);
                        sigma22 = v_sigma(si22);
                        
                        k1 = -sigma11*vk1(i1) + sigma12*vk1(i2);
                        k2 = -sigma21*vk2(i1) + sigma22*vk2(i2);
                        
                        my_sum = my_sum + ...
                            A(r2,i2,si12,si22)*...
                            conj(A(r1,i1,si11,si21))*...
                            exp(1j*(k1*x1+k2*x2));
                        
                    end
                end
            end
        end
    end
end

res = my_sum;

end

function res = fn_int(momentum)
L = 1;
k=[1.0*pi 5.0/3*pi];
k10=k(1,1);k20=k(1,2);

k11=(k10+3*k20)/2; k21=(k10-k20)/2; k12=(k10-3*k20)/2; k22=(k10+k20)/2;
%A??????????A0+??A1+??A2+??A0-??A1-??A2-
A=[1,1,0,0,0,0;0,0,0,-exp(1i*(k20-k10)*L),-exp(1i*(k21-k11)*L),0;-exp(-1i*k20*L),...
    0,1,0,0,0;0,0,0,exp(-1i*k10*L),0,-exp(1i*(k22-k12)*L);0,-exp(-1i*k21*L),-exp(1i*k12*L),...
    0,0,0;0,0,0,0,exp(-1i*k11*L),exp(1i*k22*L)];
B1=[A;1,0,0,-1,0,0];
bs=null(B1);
%A1j1=Aj+(+,+) A1j2=Aj+(+,-) A1j3=Aj+(-,+) A1j4=Aj+(-,-) A2j1=Aj-(+,+) A2j2=Aj-(+,-) A2j3=Aj-(-,+) A2j4=Aj-(-,-)
A101=bs(1);A111=bs(2);A121=bs(3);A201=bs(4);A211=bs(5);A221=bs(6);
A102=-exp(-1i*k20*L)*bs(1);A112=-exp(-1i*k21*L)*bs(2);A122=-exp(-1i*k22*L)*bs(3);
A202=-exp(1i*k20*L)*bs(4);A212=-exp(1i*k21*L)*bs(5);A222=-exp(1i*k22*L)*bs(6);
A103=-exp(1i*k10*L)*bs(1);A113=-exp(1i*k11*L)*bs(2);A123=-exp(1i*k12*L)*bs(3);
A203=-exp(-1i*k10*L)*bs(4);A213=-exp(-1i*k11*L)*bs(5);A223=-exp(-1i*k12*L)*bs(6);
A104=exp(1i*k10*L-1i*k20*L)*bs(1);A114=exp(1i*k11*L-1i*k21*L)*bs(2);
A124=exp(1i*k12*L-1i*k22*L)*bs(3);A204=exp(1i*k20*L-1i*k10*L)*bs(4);
A214=exp(1i*k21*L-1i*k11*L)*bs(5);A224=exp(1i*k22*L-1i*k12*L)*bs(6);
% [x1,x2] = meshgrid(-0.5:.02:0.5, -0.5:.02:0.5);
% f=(A101*exp(i*k10*x1+i*k20*x2)+A102*exp(i*k10*x1-i*k20*x2)+A103*exp(-i*k10*x1+i*k20*x2)+A104*exp(-i*k10*x1-i*k20*x2)
% +A111*exp(i*k11*x1+i*k21*x2)+A112*exp(i*k11*x1-i*k21*x2)+A113*exp(-i*k11*x1+i*k21*x2)+A114*exp(-i*k11*x1-i*k21*x2)
% +A121*exp(i*k12*x1+i*k22*x2)+A122*exp(i*k12*x1-i*k22*x2)+A123*exp(-i*k12*x1+i*k22*x2)+A124*exp(-i*k12*x1-i*k22*x2)).*(x2-x1<0)
% +(A201*exp(i*k10*x1+i*k20*x2)+A202*exp(i*k10*x1-i*k20*x2)+A203*exp(-i*k10*x1+i*k20*x2)+A204*exp(-i*k10*x1-i*k20*x2)
% +A211*exp(i*k11*x1+i*k21*x2)+A212*exp(i*k11*x1-i*k21*x2)+A213*exp(-i*k11*x1+i*k21*x2)+A214*exp(-i*k11*x1-i*k21*x2)
% +A221*exp(i*k12*x1+i*k22*x2)+A222*exp(i*k12*x1-i*k22*x2)+A223*exp(-i*k12*x1+i*k22*x2)+A224*exp(-i*k12*x1-i*k22*x2)).*(x2-x1>=0)

Num = 501;
vx1 = linspace(-0.5,0.5,Num);
vx2 = linspace(-0.5,0.5,Num);
[x1,x2] = meshgrid(vx1,vx2);

delta_x = 1/(Num-1);
m_nt = ones(Num,Num);
m_nt(1,1) = 0.25;
m_nt(1,Num) = 0.25;
m_nt(Num,1) = 0.25;
m_nt(Num,Num) = 0.25;
m_nt(1,2:(Num-1)) = 0.5;
m_nt(Num,2:(Num-1)) = 0.5;
m_nt(2:(Num-1),1) = 0.5;
m_nt(2:(Num-1),Num) = 0.5;

m_nt = m_nt * delta_x^2;

f=(A101*exp(1i*k10*x1+1i*k20*x2)+A102*exp(1i*k10*x1-1i*k20*x2)+A103*exp(-1i*k10*x1+1i*k20*x2)+A104*exp(-1i*k10*x1-1i*k20*x2)+...
    A111*exp(1i*k11*x1+1i*k21*x2)+A112*exp(1i*k11*x1-1i*k21*x2)+A113*exp(-1i*k11*x1+1i*k21*x2)+A114*exp(-1i*k11*x1-1i*k21*x2)+...
    A121*exp(1i*k12*x1+1i*k22*x2)+A122*exp(1i*k12*x1-1i*k22*x2)+A123*exp(-1i*k12*x1+1i*k22*x2)+A124*exp(-1i*k12*x1-1i*k22*x2)).*(x2-x1<0)+...
    (A201*exp(1i*k10*x1+1i*k20*x2)+A202*exp(1i*k10*x1-1i*k20*x2)+A203*exp(-1i*k10*x1+1i*k20*x2)+A204*exp(-1i*k10*x1-1i*k20*x2)+...
    A211*exp(1i*k11*x1+1i*k21*x2)+A212*exp(1i*k11*x1-1i*k21*x2)+A213*exp(-1i*k11*x1+1i*k21*x2)+A214*exp(-1i*k11*x1-1i*k21*x2)+...
    A221*exp(1i*k12*x1+1i*k22*x2)+A222*exp(1i*k12*x1-1i*k22*x2)+A223*exp(-1i*k12*x1+1i*k22*x2)+A224*exp(-1i*k12*x1-1i*k22*x2)).*(x2-x1>=0);

% F=abs(f);
% F1=real(f);
res = sum(sum(conj(f).*f.*exp(-1j*momentum*x1).*m_nt));

end
% surf(x1,x2,F1)
% gca=pcolor(x1,x2,F1);
% set(gca,'lineStyle','none')
% xlim([-0.5 0.5]);
% ylim([-0.5 0.5]);
% axis equal


%  subplot(4,4,1);title('c=-1');
%  subplot(4,4,2);title('c=0.1');
%  subplot(4,4,3);title('c=1');
%   subplot(4,4,4);title('c=10');
%  subplot(4,4,5);
%   subplot(4,4,6);
%  subplot(4,4,7);
%  subplot(4,4,8);
%   subplot(4,4,9);
%  subplot(4,4,10);
%    subplot(4,4,11);
%  subplot(4,4,12);
%  subplot(4,4,13);
%   subplot(4,4,14);
%  subplot(4,4,15);
%  subplot(4,4,16);
%  subplot(3,3,4);title('c=3');
%  subplot(3,3,5);title('c=4');
%  subplot(3,3,6);title('c=6');
%  subplot(3,3,7);title('c=8');
%  subplot(3,3,8);title('c=10');
%  subplot(3,3,9);title('c=100');

% theta=0:0.01:2*pi;
% Circle=zeros(2,length(theta));
% x=0;y=0;r=sqrt(3.21431544775981^2+1.79488885575112^2);
% Circle(1,:)=x+r*cos(theta);
% Circle(2,:)=y+r*sin(theta);
% Circle=Circle'



function res = fn_g(vk)

if length(vk)>1
    [kmax,max_ind] = max(vk);
    [kmin,min_ind] = min(vk);
    
    if kmax-kmin<1e-6
        res = exp(1j*mean(vk));
    else
        vk1 = vk; vk2 = vk;
        vk1(min_ind) = []; vk2(max_ind) = [];
        res = (fn_g(vk1)-fn_g(vk2))/1j/(kmax-kmin);
    end
    
else
    res = exp(1j*vk(1));
end

end


