figure(1)
hold on;
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues =0:0.1:4;
grid on;
ax.XMinorGrid = 'on';
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues =0:0.1:4;
grid on;
ax.YMinorGrid = 'on';