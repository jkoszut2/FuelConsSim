% Plot Points
figure
p1 = plot(factor1,FSAELpts(:,1),'b--');
hold on
p2 = plot(factor1,FSAELpts(:,2),'m--');
hold on
p3 = plot(factor1,FSAELpts(:,1)+FSAELpts(:,2),'r--', 'linewidth', 2);
hold on
p4 = plot(factor1,FSAEMpts(:,1),'b');
hold on
p5 = plot(factor1,FSAEMpts(:,2),'m');
hold on
p6 = plot(factor1,FSAEMpts(:,1)+FSAEMpts(:,2), 'r', 'linewidth', 2);
FSAELopt = find(FSAELpts(:,1)+FSAELpts(:,2)==max(FSAELpts(:,1)+FSAELpts(:,2)));
FSAEMopt = find(FSAEMpts(:,1)+FSAEMpts(:,2)==max(FSAEMpts(:,1)+FSAEMpts(:,2)));
p7 = xline(factor1(FSAELopt),'--');
p8 = xline(factor1(FSAEMopt));
legend([p1 p2 p3 p7 p4 p5 p6 p8], ...
       {'FSAEL Endurance', 'FSAEL Efficiency', 'FSAEL Net', 'FSAEL Maximum', ...
        'FSAEM Endurance', 'FSAEM Efficiency', 'FSAEM Net', 'FSAEM Maximum'}, ...
        'location', 'best')
title('Parameter Sweep')
xlabel('factor1')
ylabel('Points Gain')
grid on
ax1 = gca;
% Avoid scientfic notation on axes
% ax = gca;
% ax.XRuler.Exponent = 0;
%% Plot Times and Fuel Consumption
figure
p9 = plot(factor1,FSAELresults(:,1)-FSAELresults(1,1),'b--');
yyaxis right
p10 = plot(factor1,FSAELresults(:,2)-FSAELresults(1,2),'m--');
yyaxis left
hold on
p11 = plot(factor1,FSAEMresults(:,1)-FSAEMresults(1,1),'b');
yyaxis right
hold on
p12 = plot(factor1,FSAEMresults(:,2)-FSAEMresults(1,2),'m');
p13 = xline(factor1(FSAELopt),'--');
p14 = xline(factor1(FSAEMopt));
legend([p9 p10 p13 p11 p12 p14], ...
       {'FSAEL Time [sec]', 'FSAEL Fuel [cc]', 'FSAEL Maximum',...
        'FSAEM Time [sec]', 'FSAEM Fuel [cc]', 'FSAEM Maximum'}, ...
        'location', 'best')
title('Parameter Sweep')
xlabel('factor1')
ylabel('Overall Fuel Consumed Delta [cc]','rotation',-90,'VerticalAlignment','bottom')
yyaxis left
ylabel('Overall Time Delta [sec]')
grid on
ax2 = gca;
ax2.YAxis(2).Color = 'k';
% Avoid scientfic notation on axes
% ax.XRuler.Exponent = 0;