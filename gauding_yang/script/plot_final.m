close(figure(10));
figure(10);
hold on;

gca = pcolor(m_c,m_T0,m_r);
set(gca,'linestyle','none');
contour(m_c,m_T0,m_r,20,'color','yellow','linestyle','--');