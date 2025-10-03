function saveFiguresAsFrames(outputDir, figHandles, ks, varargin)
%SAVEFIGURESASFRAMES 将一组 figure 保存为有序 PNG 帧图像
%
% 用法：
%   saveFiguresAsFrames('my_frames', [fig1, fig2, fig3])
%
% 可选参数（Name-Value 对）：
%   'XLim'        - x 轴范围，如 [0 10]（默认：[]，保留原图）
%   'YLim'        - y 轴范围，如 [0 5]（默认：[]）
%   'Resolution'  - 输出分辨率（DPI，默认 150）
%   'BackgroundColor' - 背景色（默认 'white'）
%   'Prefix'      - 文件名前缀（默认 'frame'）
%
% 输出文件：outputDir/frame001.png, frame002.png, ...

    % === 默认参数 ===
    p = inputParser;
    addRequired(p, 'outputDir', @ischar);
    addRequired(p, 'figHandles', @(x) isvector(x) && all(isvalid(x)) && all(strcmp({x.Type}, 'figure')));
    addParameter(p, 'XLim', [0 1], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
    addParameter(p, 'YLim', [0 1], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
    addParameter(p, 'Resolution', 150, @(x) isnumeric(x) && x > 0);
    addParameter(p, 'BackgroundColor', 'white', @ischar);
    addParameter(p, 'Prefix', 'frame', @ischar);
    parse(p, outputDir, figHandles, varargin{:});

    xlim_range = p.Results.XLim;
    ylim_range = p.Results.YLim;
    resolution = p.Results.Resolution;
    bg_color   = p.Results.BackgroundColor;
    prefix     = p.Results.Prefix;

    % === 创建输出目录 ===
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    % === 遍历每个 figure ===
    nFigs = numel(figHandles);
    for k = 1:nFigs
        fig = figHandles(k);
        
        % 设置背景色
        fig.Color = bg_color;
        
        % 获取第一个 axes（如有）
        axList = findobj(fig, 'Type', 'axes');
        if ~isempty(axList)
            ax = axList(1); % 取第一个坐标轴
            if ~isempty(xlim_range)
                xlim(ax, xlim_range);
            end
            if ~isempty(ylim_range)
                ylim(ax, ylim_range);
            end
        end

        % 构造文件名
        filename = fullfile(outputDir, sprintf('%s%03d.png', prefix, k + ks));

        % 保存图像
        if verLessThan('matlab', '9.8') % R2020a 之前
            print(fig, '-dpng', ['-r' num2str(resolution)], filename);
        else
            exportgraphics(fig, filename, ...
                'Resolution', resolution, ...
                'ContentType', 'auto', ...
                'BackgroundColor', bg_color);
        end

        fprintf('Saved: %s\n', filename);
    end

    fprintf('✅ 共 %d 个 figure 已保存到 "%s"\n', nFigs, outputDir);
end