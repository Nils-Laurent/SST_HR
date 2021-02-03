function save_fig2_EMD_full(T_hsz, BPM_X, BPM_comp,...
    EMD_STFT, EMD_SST, name)

axisFSZ = 22;
labelSZ = 36;
lgdSZ = 28;
lenSZ = 700;

figure;
imagesc(T_hsz, BPM_comp, EMD_STFT);
axis xy;
colormap(gray);
ax = gca;
ax.XAxis.FontSize = axisFSZ;
ax.YAxis.FontSize = axisFSZ;
xlabel('time', 'interpreter', 'latex', 'FontSize', labelSZ);
ylabel('component bpm', 'interpreter', 'latex', 'FontSize', labelSZ);
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, lenSZ, lenSZ])
saveas(gcf, strcat("fig2_EMD_full_STFT_", name),'epsc');

figure;
imagesc(T_hsz, BPM_comp, EMD_SST);
axis xy;
colormap(gray);
ax = gca;
ax.XAxis.FontSize = axisFSZ;
ax.YAxis.FontSize = axisFSZ;
xlabel('time', 'interpreter', 'latex', 'FontSize', labelSZ);
ylabel('component bpm', 'interpreter', 'latex', 'FontSize', labelSZ);
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, lenSZ, lenSZ])
saveas(gcf, strcat("fig2_EMD_full_SST_", name),'epsc');

end

