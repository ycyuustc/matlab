
v_h_07 = v_h;
v_t_07 = v_t;
v_mu_07 = v_mu;
save finished_06.mat
% v_h_05 = v_h;
% v_t_05 = v_t;
% v_mu_05 = v_mu;
% save finished_05.mat

% v_h_04 = v_h;
% v_t_04 = v_t;
% v_mu_04 = v_mu;
% save finished_04.mat
% v_h_02 = v_h;
% v_t_02 = v_t;
% v_mu_02 = v_mu;

% v_h_01 = v_h;
% v_t_01 = v_t;
% v_mu_01 = v_mu;

% v_h_03 = v_h;
% v_t_03 = v_t;
% v_mu_03 = v_mu;

% v_h_04 = v_h;
% v_t_04 = v_t;
% v_mu_04 = v_mu;

% mu0 = -0.2427852; % for t0 = 0.002
% mu0 = -0.242589;  % for t0 = 0.003
%  mu0 = -0.242347; % for t0 = 0.004
% mu0 = -0.2421007; % for t0 = 0.005
% mu0 = -0.2418725; % for 0.006
%  mu0 = -0.241672; % 0.007
mu0 = -0.241499; % 0.008
h0 = 0.3;
t0 = 0.008;

tv = pl_gaudin_insert_epsilon(h0,mu0,t0);
n0 = tv(2);
s0 = tv(6);

v_h1 = linspace(0.3,0.47,20);
v_h2 = linspace(0.47,0.51,50);
v_h3 = linspace(0.51,0.61,20);
v_h4 = linspace(0.61,0.70,50);
v_h5 = linspace(0.70,0.9,20);
v_h = [v_h1(1:end-1),v_h2(1:end-1),v_h3(1:end-1),v_h4(1:end-1),v_h5];
Num = length(v_h);
disp(Num);

v_mu = zeros(1,Num);
v_t = zeros(1,Num);

mu = mu0;
t = t0;
for ih = 1:Num
    
   h = v_h(ih);
   disp(fn_judge(h,mu,t,n0,s0));
   [mu,t] = fn_find(h,mu,t,n0,s0);
   disp(fn_judge(h,mu,t,n0,s0));
   v_mu(ih) = mu;
   v_t(ih) = t;
   disp(ih);
   disp('=================================');
   
end




%%%%%%%%%%%%%%%%%%%%% functions %%%%%%%%%%%%%%%%%%%%%%

function [mu,t] = fn_find(h0,mu0,t0,n0,s0)

mu = mu0;
t = t0;
delta = 0.01;

while delta>1e-7
    
   v_mu = [mu-delta,mu,mu+delta];
   v_t = [t-delta*0.1,t,t+delta*0.1];
   m_r = zeros(3,3);
   for i_mu = 1:3
       for i_t = 1:3
            tmu = v_mu(i_mu);
            tt = v_t(i_t);
            m_r(i_mu,i_t) = fn_judge(h0,tmu,tt,n0,s0);
       end
   end
   
   [ind_mu,ind_t] = find(m_r==min(min(m_r)));
   mu = v_mu(ind_mu);
   t = v_t(ind_t);
   if ind_mu==2 && ind_t==2
       delta = delta/2;
       disp(delta);
   end
    
end

end



function result = fn_judge(h,mu,t,n0,s0)

tv = pl_gaudin_insert_epsilon(h,mu,t);
% n = pl_gaudin_density(h,mu,t);
% s = pl_gaudin_entropy(h,mu,t);
n = tv(2);
s = tv(6);

result = sqrt((n-n0)^2+(s-s0)^2);

end