th_str = 'Posteriori for theta when Psi = '; th_s = num2str(psi);
l_str = 'Posteriori for lambda when Psi = ';
t_str = 'Posteriori for the t when Psi = '; 
figure
subplot(3,1,1)
hist(th_long,50)
ylim([0 800])
title([th_str, th_s])
subplot(3,1,2)
hist(l_long, 50)%,round(max(l_long)*8))
%xlim([0 5])
%ylim([0 7000])
title([l_str, th_s])
subplot(3,1,3)
hist(t_long)
title([t_str, th_s])
ylim([0 12000])

% max(l_long(:))