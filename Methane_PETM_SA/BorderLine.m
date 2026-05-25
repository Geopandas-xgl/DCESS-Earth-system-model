function BorderLine(ax, opts)
    arguments
        ax matlab.graphics.axis.Axes
        opts.LineWidth (1,1) double = 1
        opts.Side (1,1) string {mustBeMember(opts.Side,["","xre","yre"])} = ""
    end

    % 获取极限
    xL = ax.XLim; yL = ax.YLim;
    hold(ax,'on');

    % 默认画右侧和顶部边线
    if opts.Side == ""
        line(ax, [xL(2) xL(2)], yL, 'Color','k','LineWidth',opts.LineWidth);
        line(ax, xL, [yL(2) yL(2)], 'Color','k','LineWidth',opts.LineWidth);
    end

    % left side
    if opts.Side == "xre"
        line(ax, [xL(1) xL(1)], yL, 'Color','k','LineWidth',opts.LineWidth);
        line(ax, xL, [yL(2) yL(2)], 'Color','k','LineWidth',opts.LineWidth);
    end

    % top side
    if opts.Side == "yre"
        line(ax, xL, [yL(1) yL(1)], 'Color','k','LineWidth',opts.LineWidth);
        line(ax, [xL(2) xL(2)], yL, 'Color','k','LineWidth',opts.LineWidth);
    end

    hold(ax,'off');
end
