function legend_Ismall(loc)


lgdSZ = 24;
lgd = legend;
lgd.FontSize = lgdSZ;
if nargin > 0
    set(lgd,'Interpreter','latex', 'Location', loc);
else
    set(lgd,'Interpreter','latex');
end
end

