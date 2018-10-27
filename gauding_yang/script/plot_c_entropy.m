v_c = linspace(0.8, 1.2, 1000);
v_entropy = zeros(1,1000);

h0 = 0.5;
mu0 = -0.2;
T0 = 0.0001;

for k=1:1000
   
    coupling = v_c(k);
    h = h0/coupling^2;
    mu = mu0/coupling^2;
    T = T0/coupling^2;
    
    v_entropy(k) = pl_gaudin_entropy(h,mu,T)*coupling;
    
    disp(k);
    
end



