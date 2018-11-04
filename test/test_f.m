function res = test_f(x)
global x1;
global x2;
global A0;
global A1;
global A2;
res=((x-x1-1)*(x-x2-1)/(x-x1+1)/(x-x2+1)/A0-1)^2 ...
+((x1-x-1)/(x1-x+1)/A1-1)^2 ...
+((x2-x-1)/(x2-x+1)/A2-1)^2;
end