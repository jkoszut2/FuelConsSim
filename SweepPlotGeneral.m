% factor1 = factor1-13500;
plot(factor1,FSAELpts(:,1),'b--')
hold on
plot(factor1,FSAELpts(:,2),'m--')
hold on
plot(factor1,FSAELpts(:,1)+FSAELpts(:,2),'r--', 'linewidth', 2)
hold on
plot(factor1,FSAEMpts(:,1),'b')
hold on
plot(factor1,FSAEMpts(:,2),'m')
hold on
plot(factor1,FSAEMpts(:,1)+FSAEMpts(:,2), 'r', 'linewidth', 2)
legend('FSAEL end', 'FSAEL eff', 'FSAEL Net', ...
       'FSAEM end', 'FSAEM eff', 'FSAEM Net', 'location', 'best')
% title('(FEPW and Torque) * factor1')
title('Parameter Sweep')
xlabel('factor1')
ylabel('pts gain')
grid on
% Avoid scientfic notation on axes
% ax = gca;
% ax.XRuler.Exponent = 0;