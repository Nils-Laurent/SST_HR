function save_fig3_std(SNR_in, std_vec_SST, std_vec_STFT, name)

axisFSZ = 20;
labelSZ = 32;
lgdSZ = 24;
lenSZ = 700;

figure;
hold on;
plot(SNR_in, std_vec_SST, 'linewidth', 2, 'DisplayName', 'std SST');
plot(SNR_in, std_vec_STFT, '--', 'linewidth', 2, 'DisplayName', 'std STFT');
hold off;
lgd = legend;
lgd.FontSize = lgdSZ;
xlim([SNR_in(1), SNR_in(end)]);

xlabel('SNR in', 'interpreter', 'latex', 'FontSize', labelSZ);
ylabel('Standard deviation', 'interpreter', 'latex', 'FontSize', labelSZ);
ax = gca;
ax.XAxis.FontSize = axisFSZ;
ax.YAxis.FontSize = axisFSZ;
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, lenSZ, lenSZ])
saveas(gcf, strcat("fig3_std_", name),'epsc');

end

