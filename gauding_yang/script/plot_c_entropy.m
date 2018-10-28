close(figure(1));
figure(1);
hold on;

v_c = linspace(0.8, 1.2, 1000);
v_entropy = zeros(1,1000);
v_density = zeros(1,1000);
v_mag = zeros(1,1000);

h0 = 0.5;
mu0 = -0.2;
T0 = 0.001;

for k=1:1000
   
    coupling = v_c(k);
    h = h0/coupling^2;
    mu = mu0/coupling^2;
    T = T0/coupling^2;
    
    v_entropy(k) = pl_gaudin_entropy(h,mu,T)*coupling;
%     v_density(k) = pl_gaudin_density(h,mu,T)*coupling;
%     v_mag(k) = pl_gaudin_mag(h,mu,T)*coupling;
    disp(k);
    
end

plot(v_c,v_entropy,'-*');
legend('T=1e-3');




T0 = 0.002;
for k=1:1000
   
    coupling = v_c(k);
    h = h0/coupling^2;
    mu = mu0/coupling^2;
    T = T0/coupling^2;
    
    v_entropy(k) = pl_gaudin_entropy(h,mu,T)*coupling;
    disp(k);
    
end

plot(v_c,v_entropy,'-*');
legend('T=2e-3');



T0 = 0.003;
for k=1:1000
   
    coupling = v_c(k);
    h = h0/coupling^2;
    mu = mu0/coupling^2;
    T = T0/coupling^2;
    
    v_entropy(k) = pl_gaudin_entropy(h,mu,T)*coupling;
    disp(k);
    
end

plot(v_c,v_entropy,'-*');
legend('T=3e-3');


T0 = 0.004;
for k=1:1000
   
    coupling = v_c(k);
    h = h0/coupling^2;
    mu = mu0/coupling^2;
    T = T0/coupling^2;
    
    v_entropy(k) = pl_gaudin_entropy(h,mu,T)*coupling;
    disp(k);
    
end

plot(v_c,v_entropy,'-*');
legend('T=4e-3');