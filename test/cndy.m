% X = @(k1,k2,g) (k2-k1+1j*g).*(k2+k1-1j*g)./(k2-k1-1j*g)./(k2+k1+1j*g);
% f1 = @(k1,k2,g,L,d) (1+X(k1,k2,g).*exp(1j*k1*L))./(1-X(k1,k2,g).*exp(1j*k1*L))...
%     -(exp(1j*k1*L)+1)./(exp(1j*k1*L)-1)-2*d./1j./k1;
% 
% f2 = @(k1,k2,g,L,d) 1./tan(k1*L/2+atan(2*k1*g./(g^2+k2^2-k1^2)))+1./tan(k1*L/2)+2*d./k1;
% 
% f1(1,2,3,4,5)
% f2(1,2,3,4,5)
% 
% vk1 = linspace(-2*pi,2*pi,1001);
% vk2 = linspace(-2*pi,2*pi,1001);
% [mk1,mk2] = meshgrid(vk1,vk2);
% 
% 
% theta = @(k1,k2,g) atan(2*k1*g./(g^2+k2.^2-k1.^2));
% 
% f3 = @(k1,k2,g,L,d) k1.*sin(k1*L+theta(k1,k2,g))...
%     + 2*d.*sin(k1*L/2).*sin(k1*L/2+theta(k1,k2,g));
% figure(1);
% contour(mk1,mk2,(f3(mk1,mk2,1,2,0.5)),[0.00001,0.1]);
% hold on;
% contour(mk1,mk2,(f3(mk2,mk1,1,2,0.5)),[0.00001,0.1]);
% 
% 
% 
% 
% f4 = @(k1,k2,g,L,d) sin(k1*L/2+theta(k1,k2,g)).^2.*sin(k2*L/2).^2 ...
%     -sin(k2*L/2+theta(k2,k1,g)).^2.*sin(k1*L/2).^2;
% figure(1);
% contour(mk1,mk2,(f4(mk1,mk2,1,2,0.5)),[0.0000001,0.1]);
% 
f5 = @(k1,k2,g,L,d) (g^2+k2.^2-k1.^2).*(k1.*sin(k1*L)+d.*(1-cos(k1*L))) ...
    +2.*k1.*g.*(k1.*cos(k1*L)+d*sin(k1*L));
% 
% figure(5);
% contour(mk1,mk2,(f5(mk1,mk2,1,2,0.5)),[-0.000001,0.001]);
% hold on;
% contour(mk1,mk2,(f5(mk2,mk1,1,2,0.5)),[-0.000001,0.001]);
% contour(mk1,mk2,(f4(mk1,mk2,1,2,0.5)),[0.0000001,0.001]);
% 
% figure(10);
% figure(10);hold on;
% contour(mk1,mk2,(f4(mk1,mk2,1,2,0.5)),[0.0000001,0.001]);
% 
% 
% figure(6);
% contour(mk1,mk2,(f3(mk1,mk2,1,2,0.5)),[0.00001,0.1]);
% 
% th = @(k1,k2,g,L) k1*L/2+atan(2*k1*g./(g^2+k2.^2-k1.^2));
% 
% f8 = @(k1,k2,g,L,d) k1.*sin(th(k1,k2,g,L)+k1*L/2).*sin(th(k2,k1,g,L)).*sin(k2*L/2)...
%     -k2.*sin(th(k2,k1,g,L)+k2*L/2).*sin(th(k1,k2,g,L)).*sin(k1*L/2);
% 
% figure(8);
% contour(mk1,mk2,f8(mk1,mk2,1,2,0.5),[0.000001,0.0001]);
% 
% f9 = @(k1,k2,g,L,d) sin(th(k1,k2,g,L)).*sin(k2*L/2) ...
%     +sin(th(k2,k1,g,L)).*sin(k1*L/2);
% figure(9);
% hold on;
% contour(mk1,mk2,f9(mk1,mk2,1,2,0.5),[0.000001,0.0001]);
% plot(k1,k2,'r*');

f10 = @(k1,k2,g,L,d) k1.^2.*(d*(1-cos(k2*L))+k2.*sin(k2*L)) ...
    -k2.^2.*(d*(1-cos(k1*L))+k1.*sin(k1*L));

% vk1 = linspace(kk1-0.00001,kk1+0.00001,1000);
% vk2 = linspace(kk2-0.00001,kk2+0.00001,1000);
% 
% [mk1,mk2] = meshgrid(vk1,vk2);
% 
% tm = abs(f5(mk1,mk2,g,L,d))+abs(f5(mk2,mk1,g,L,d));

vk1 = linspace(-8*pi,8*pi,1001);
vk2 = linspace(-4*pi,4*pi,1001);
[mk1,mk2] = meshgrid(vk1,vk2);

g=1.3;L=1;d=0.58;
close(figure(1));
figure(1);
contour(mk1,mk2,f5(mk1,mk2,g,L,d),[0.000001,0.0001]);
hold on;
contour(mk1,mk2,f5(mk2,mk1,g,L,d),[0.000001,0.0001]);
contour(mk1,mk2,f10(mk1,mk2,g,L,d),[0.000001,0.0001]);

% vx = linspace(0.01,10*pi,1000);
% figure(888);
% plot(vx,(d*(1-cos(vx*L)+vx.*sin(vx*L)))./vx.^2,'r*');



