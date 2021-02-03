function write_fig_TF(xV, label_x, yV, label_y, Im, name)

if isempty(xV)
    xV = 1:size(Im, 2);
end
if isempty(yV)
    yV = 1:size(Im, 1);
end

figure;
imagesc(xV, yV, Im);
axis xy;

xlabel(label_x, 'interpreter', 'latex');
ylabel(label_y, 'interpreter', 'latex');
colormap(flipud(gray));
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, 1000, 1000])

savefig(name);
saveas(gcf, name, 'epsc');
close all;

end