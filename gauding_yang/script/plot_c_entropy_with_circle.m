close(figure(1));
figure(1);
h1 = pcolor(m_c,m_T0,m_r);
set(h1,'linestyle','none');

th2x = h2x(h2x>0.85);
th2y = h2y(h2x>0.85);
th8x = h8x(h8x>0.85);
th8y = h8y(h8x>0.85);

hold on;
plot(th2x,th2y,'--w');
plot(th8x,th8y,'--w');
plot([th2x(1),th8x(1)],[th2y(1),th8y(1)],'--w');
plot([th2x(end),th8x(end)],[th2y(end),th8y(end)],'--w');