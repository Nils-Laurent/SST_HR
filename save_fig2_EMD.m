function save_fig2_EMD(T_hsz, BPM_X, BPM_comp,...
    R_STFT, R_SST, name)

axisFSZ = 20;
labelSZ = 32;
lgdSZ = 24;
lenSZ = 700;

figure;
imagesc(T_hsz, BPM_X, R_STFT.ASTFT);
axis xy;
colormap(flipud(gray));
ax = gca;
ax.XAxis.FontSize = axisFSZ;
ax.YAxis.FontSize = axisFSZ;
xlabel('time', 'interpreter', 'latex', 'FontSize', labelSZ);
ylabel('bpm', 'interpreter', 'latex', 'FontSize', labelSZ);
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, lenSZ, lenSZ])
saveas(gcf, strcat("fig2_STFT_", name),'epsc');

figure;
imagesc(T_hsz, BPM_X, R_SST.ASST);
axis xy;
colormap(flipud(gray));
ax = gca;
ax.XAxis.FontSize = axisFSZ;
ax.YAxis.FontSize = axisFSZ;
xlabel('time', 'interpreter', 'latex', 'FontSize', labelSZ);
ylabel('bpm', 'interpreter', 'latex', 'FontSize', labelSZ);
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, lenSZ, lenSZ])
saveas(gcf, strcat("fig2_SST_", name),'epsc');

figure;
imagesc(T_hsz, BPM_comp, R_STFT.EMD);
axis xy;
colormap(gray);
hold on;
plot(T_hsz, BPM_comp(R_STFT.LB), '--', 'linewidth', 2, 'DisplayName', 'Lower bound');
plot(T_hsz, BPM_comp(R_STFT.HB), '-', 'linewidth', 2, 'DisplayName', 'Upper bound');
plot(T_hsz, BPM_comp(R_STFT.CVec), 'g-', 'DisplayName', 'HR detection');
hold off;

lgd = legend;
lgd.FontSize = lgdSZ;

ax = gca;
ax.XAxis.FontSize = axisFSZ;
ax.YAxis.FontSize = axisFSZ;
xlabel('time', 'interpreter', 'latex', 'FontSize', labelSZ);
ylabel('component bpm', 'interpreter', 'latex', 'FontSize', labelSZ);
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, lenSZ, lenSZ])
saveas(gcf, strcat("fig2_EMD_STFT_", name),'epsc');

figure;
imagesc(T_hsz, BPM_comp, R_SST.EMD);
axis xy;
colormap(gray);
hold on;
plot(T_hsz, BPM_comp(R_SST.LB), '--', 'linewidth', 2, 'DisplayName', 'Lower bound');
plot(T_hsz, BPM_comp(R_SST.HB), '-', 'linewidth', 2, 'DisplayName', 'Upper bound');
plot(T_hsz, BPM_comp(R_SST.CVec), 'g-', 'DisplayName', 'HR detection');
hold off;

lgd = legend;
lgd.FontSize = lgdSZ;

ax = gca;
ax.XAxis.FontSize = axisFSZ;
ax.YAxis.FontSize = axisFSZ;
xlabel('time', 'interpreter', 'latex', 'FontSize', labelSZ);
ylabel('component bpm', 'interpreter', 'latex', 'FontSize', labelSZ);
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, lenSZ, lenSZ])
saveas(gcf, strcat("fig2_EMD_SST_", name),'epsc');

end

