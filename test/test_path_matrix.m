a = 3;
b = 1;

Num = 8;

m = diag(a*ones(1,Num))+diag(b*ones(1,Num-1),-1)+diag(b*ones(1,Num-1),1);

v_eig = a + b*2*cos(2*pi/(2*Num+2)*(1:Num));
