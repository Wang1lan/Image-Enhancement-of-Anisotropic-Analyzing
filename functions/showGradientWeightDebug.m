function figs = showGradientWeightDebug(debug, displayOpts)
%SHOWGRADIENTWEIGHTDEBUG 展示 maps 级梯度分析，默认不创建目录或保存。
%   figs = showGradientWeightDebug(debug,struct('save',false))
%   figs.maps 为图窗；保存时 figs.outputTag 返回该组导出的公共文件名前缀。
%   save=true 时指定 outputDir/prefix，输出分析 PNG、G/W 统计 CSV、精确 G/W 参数 MAT。
%   G 图裁剪色阶仅影响显示；W 图固定为 [0,1]，计算不使用显示归一化。

required = {'Yinput','Yguide','G','W','validMask','paramsUsed','stats'};
assert(isstruct(debug) && isscalar(debug) && all(isfield(debug, required)), ...
    'showGradientWeightDebug:MapsRequired', '请传入 computeGradientWeight 的 maps 级 debug。');
if nargin < 2
    displayOpts = struct();
end
if ~isfield(displayOpts, 'save')
    displayOpts.save = false;
end
if displayOpts.save
    assert(all(isfield(displayOpts, {'outputDir','prefix'})) && ...
        ~isempty(displayOpts.outputDir) && ~isempty(displayOpts.prefix), ...
        'showGradientWeightDebug:SaveOptions', '保存时必须指定 outputDir 和 prefix。');
end
p = debug.paramsUsed;
method = 'raw';
if p.denoise.enable && p.denoise.windowSize > 1
    method = p.denoise.method;
end
label = sprintf('GW_%s | GLow=%.6g, GHigh=%.6g, Wmin=%g', ...
    method, p.GLow, p.GHigh, p.Wmin);
figs.maps = figure('Color', 'w', 'Name', ['梯度权重：' label], ...
    'NumberTitle', 'off', 'Position', [80 70 1250 780]);
layout = tiledlayout(figs.maps, 2, 3, 'TileSpacing', 'compact', 'Padding', 'loose');
maps = {debug.Yinput, debug.Yguide, debug.G, [], debug.W, double(debug.validMask)};
titles = {'分析输入 Y / inputMax', '保边引导图 Yguide', ...
    'Sobel L1 梯度 G（显示上限 GHigh）', '梯度到权重的映射', ...
    '最终权重 W', '有效区标记'};
for k = 1:6
    ax = nexttile(layout);
    if k == 4
        g = linspace(0, min(1, 1.25*p.GHigh), 256);
        w = p.Wmin + (1-p.Wmin).*min(max((g-p.GLow)/(p.GHigh-p.GLow), 0), 1);
        plot(ax, g, w, 'LineWidth', 1.5);
        hold(ax, 'on');
        xline(ax, p.GLow, '--r', 'GLow');
        xline(ax, p.GHigh, '--k', 'GHigh');
        ylim(ax, [0 1]); grid(ax, 'on');
        xlabel(ax, '归一化 Sobel L1 梯度 G'); ylabel(ax, 'W');
    else
        limits = [0 1];
        if k == 3
            limits = [0 p.GHigh];
        end
        imagesc(ax, maps{k}, limits);
        axis(ax, 'image'); axis(ax, 'off');
        if k <= 2 || k == 6
            colormap(ax, gray(256));
        else
            colormap(ax, hot(256));
        end
        colorbar(ax);
    end
    if k == 4
        title(ax, {sprintf('%s | %s | Wmin=%g', titles{k}, method, p.Wmin), ...
            sprintf('GLow=%.5g, GHigh=%.5g', p.GLow, p.GHigh)}, ...
            'Interpreter', 'none', 'FontSize', 10);
    else
        title(ax, titles{k}, 'Interpreter', 'none');
    end
end
drawnow; % 刷新不可见图窗的布局后再导出

if displayOpts.save
    outputDir = displayOpts.outputDir;
    if ~isfolder(outputDir)
        mkdir(outputDir);
    end
    tag = sprintf('%s_GW_%s_G%.6g-%.6g_W%g_%s', char(displayOpts.prefix), ...
        method, p.GLow, p.GHigh, p.Wmin, char(datetime('now', 'Format', 'yyyyMMdd''T''HHmmssSSS')));
    figs.outputTag = tag;
    exportgraphics(figs.maps, fullfile(outputDir, [tag '_maps.png']), 'Resolution', 150);
    records = [debug.stats.G; debug.stats.W];
    T = addvars(struct2table(records), ["G";"W"], 'Before', 1, 'NewVariableNames', 'metric');
    ratioNames = fieldnames(debug.stats.ratios);
    for k = 1:numel(ratioNames)
        name = ratioNames{k};
        T.(['ratio' name]) = repmat(debug.stats.ratios.(name), 2, 1);
    end
    writetable(T, fullfile(outputDir, [tag '_stats.csv']));
    paramsUsed = p;
    stats = debug.stats;
    G = debug.G;
    W = debug.W;
    validMask = debug.validMask;
    save(fullfile(outputDir, [tag '_params.mat']), 'paramsUsed', 'stats', 'G', 'W', 'validMask');
end
end
