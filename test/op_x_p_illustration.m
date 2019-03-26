L = 50;
N = 9;
delta_x = L/N;
delta_p = 2*pi/L;

% vx vp are bases:
vx = linspace(0,L,N+1);
vx = vx(1:end-1);
vp = linspace(-2*pi/L*(N-1)/2,2*pi/L*(N-1)/2,N);

mx = eye(N);
[m1,m2] = meshgrid(vp,vx);
mp = exp(1j*m1.*m2);
mx = mx./sqrt(delta_x);
mp = mp./sqrt(delta_p)/sqrt(N);
mxd = mx';
mpd = mp';

op_P = zeros(N,N);
op_X = zeros(N,N);
for i=1:N
    
   op_P = delta_p*vp(i)*mp(:,i)*mpd(i,:)+op_P;
   op_X = delta_x*vx(i)*mx(:,i)*mxd(i,:)+op_X;
    
end


% 
tm = zeros(N,N);
for i = 1:N
   tm = tm + delta_x*mx(:,i)*mxd(i,:); 
end
disp(max(max(abs(eye(N)-tm))));

tm = zeros(N,N);
for i = 1:N
   tm = tm + delta_p*mp(:,i)*mpd(i,:); 
end
disp(max(max(abs(eye(N)-tm))));

nx = 2;
np = 5;

tvx = mx(:,nx);
tvp = mp(:,np);

tm1 = op_P*tvx;
tm2 = zeros(N,N);
for i=1:N
   tm2 = tm2 + delta_p*vp(i)*exp(-1j*vp(i)*vx(nx))/sqrt(2*pi)*mp(:,i); 
end
disp(max(max(abs(tm1-tm2))));

tm1 = tvp;
tm2 = zeros(N,N);
for i=1:N
   tm2 = tm2 + delta_x*exp(1j*vx(i)*vp(np))/sqrt(2*pi)*mx(:,i); 
end
disp(max(max(abs(tm1-tm2))));

nx = 7;
tm1 = op_P*mx(:,nx);
tm2 = zeros(N,N);
for ix=1:N
    for ip=1:N
        
        tm2 = tm2 + delta_x*delta_p...
            *1/2/pi*exp(1j*(vx(ix)-vx(nx))*vp(ip))*vp(ip)*mx(:,ix);
        
    end
end
disp(max(max(abs(tm1-tm2))));

vx2 = linspace(-L,L,2*N+1);

[mp2,mx2] = meshgrid(vp,vx2);
tm = mp2/2/pi.*exp(1j*mx2.*mp2)*delta_p;
vdisx = sum(tm,2);

close(figure(1));figure(1);
plot(vx2,imag(vdisx),'-r*');

close(figure(2));figure(2);
plot(vx2,real(vdisx),'-b*');

vnx = -N:N;
fx = -2*pi/L*N.*exp(1j*2*pi/N*vnx*(N+1)/2)./(1-exp(1j*2*pi/N*vnx));
fx(1) = 0;
fx(N+1) = 0;
fx(2*N+1) = 0;

fx_app = 1/1j*N^2/L*power(-1,vnx)./vnx;
fx_app(1) = 0;
fx_app(N+1) = 0;
fx_app(2*N+1) = 0;

figure(1);
hold on;
plot(vx2,imag(fx),'go');
plot(vx2,imag(fx_app),'b');
figure(2);
hold on;
plot(vx2,real(fx),'go');
plot(vx2,real(fx_app),'b');

close(figure(3));
figure(3);hold on;
plot(vp,vp,'bo');
vgp = (exp(1j*delta_x*vp)-1)/1j/delta_x;
plot(vp,real(vgp),'r*');

op_Papp = zeros(N,N);
for i=1:N
    
   op_Papp = delta_p*vgp(i)*mp(:,i)*mpd(i,:)+op_Papp;
    
end

close(figure(4));
figure(4);hold on;
plot(vp,vp,'bo');
vgp2 = (exp(1j*delta_x*vp)-exp(-1j*delta_x*vp))/1j/2/delta_x;
plot(vp,real(vgp2),'r*');

op_Papp2 = zeros(N,N);
for i=1:N
    
   op_Papp2 = delta_p*vgp2(i)*mp(:,i)*mpd(i,:)+op_Papp2;
    
end

H = op_P*op_P/2;
delta_t = 1.0;

tm = expm(-1j*delta_t*H);
% disp(tm);

nx1 = 3;
nx2 = 4;
q1 = vx(nx1);q2 = vx(nx2);
disp(mxd(nx2,:)*tm*mx(:,nx1));
tv = delta_p/2/pi*exp(-1j*delta_t/2*vp.^2).*exp(1j*(q2-q1)*vp);
disp(sum(tv));

close(figure(5));figure(5);
plot(vp,real(tv),'-b*');hold on; plot(vp,imag(tv),'-r*');

tv = -delta_t*vp.^2/2+(q2-q1)*vp;
close(figure(6));figure(6);
plot(vp,tv,'-b*');

