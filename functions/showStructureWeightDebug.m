function figs = showStructureWeightDebug(debug, displayOpts)
%SHOWSTRUCTUREWEIGHTDEBUG 显示 maps 档分析图、有效区分布与门限。
%   figs = showStructureWeightDebug(debug, struct('save', false))
%   返回 figs.maps（2x3 分析图）和 figs.distributions（C/S 分布）。
%   S 始终以 log10(S+epsilon) 显示；WS 的计算仍在原始线性域。
%   保存时必须设置 save=true、outputDir、prefix；prefix 由调用方包含
%   输入图名及 ROI 坐标/标签。文件名追加 sigma/rho、四个门限和时间。
%   导出两张 PNG、统计 CSV、参数 MAT；默认不创建目录、不保存文件。

required = {'Y', 'C', 'S', 'WC', 'WS', 'W', 'validMask', 'paramsUsed', 'stats'};
if ~isstruct(debug) || ~isscalar(debug) || ~all(isfield(debug, required))
    error('showStructureWeightDebug:MapsRequired', 'Use the maps-level second output of computeStructureWeight.');
end
if nargin < 2
    displayOpts = struct();
end
if ~isfield(displayOpts, 'save')
    displayOpts.save = false;
end
if displayOpts.save && (~all(isfield(displayOpts, {'outputDir', 'prefix'})) ...
        || isempty(displayOpts.outputDir) || isempty(displayOpts.prefix))
    error('showStructureWeightDebug:SaveOptions', 'Saving requires outputDir and prefix.');
end
p = debug.paramsUsed;
figs.maps = figure('Color', 'w', 'Name', '结构权重：分析图', ...
    'NumberTitle', 'off', 'Position', [80 70 1250 780]);
maps = {debug.Y, debug.C, log10(debug.S+p.epsilon), debug.WC, debug.WS, debug.W};
titles = {'分析输入 Y', '方向一致性 C', 'log10(S + epsilon)（仅显示）', ...
    '方向权重 WC', '能量权重 WS', '最终权重 W'};
for k = 1:6
    ax = subplot(2, 3, k, 'Parent', figs.maps);
    imagesc(ax, maps{k});
    axis(ax, 'image'); axis(ax, 'off');
    if k ~= 3
        caxis(ax, [0 1]);
    end
    if k == 1
        colormap(ax, gray(256));
    else
        colormap(ax, hot(256));
    end
    colorbar(ax);
    title(ax, titles{k}, 'Interpreter', 'none');
end
sgtitle(figs.maps, sprintf('sigma=%g, rho=%g | C=[%g, %g], S=[%.5g, %.5g]', ...
    p.sigma, p.rho, p.CLow, p.CHigh, p.SLow, p.SHigh), 'Interpreter', 'none');

mask = debug.validMask;
figs.distributions = figure('Color', 'w', 'Name', '结构权重：有效区分布', ...
    'NumberTitle', 'off', 'Position', [110 110 1100 420]);
ax = subplot(1, 2, 1, 'Parent', figs.distributions);
histogram(ax, debug.C(mask), linspace(0, 1, 65), 'Normalization', 'probability');
hold(ax, 'on');
xline(ax, p.CLow, '--r', 'CLow'); xline(ax, p.CHigh, '--k', 'CHigh');
xlim(ax, [0 1]); grid(ax, 'on');
xlabel(ax, 'C'); ylabel(ax, '有效像素比例'); title(ax, '方向一致性');
ax = subplot(1, 2, 2, 'Parent', figs.distributions);
histogram(ax, log10(debug.S(mask)+p.epsilon), 64, 'Normalization', 'probability');
hold(ax, 'on');
xline(ax, log10(p.SLow+p.epsilon), '--r', 'SLow');
xline(ax, log10(p.SHigh+p.epsilon), '--k', 'SHigh');
grid(ax, 'on'); xlabel(ax, 'log10(S + epsilon)'); ylabel(ax, '有效像素比例');
title(ax, '能量分布（对数显示）');

if displayOpts.save
    outputDir = displayOpts.outputDir;
    if ~isfolder(outputDir)
        mkdir(outputDir);
    end
    tag = sprintf('%s_s%g_r%g_C%g-%g_S%g-%g_%s', char(displayOpts.prefix), ...
        p.sigma, p.rho, p.CLow, p.CHigh, p.SLow, p.SHigh, datestr(now, 'yyyymmddTHHMMSSFFF'));
    exportgraphics(figs.maps, fullfile(outputDir, [tag '_maps.png']), 'Resolution', 150);
    exportgraphics(figs.distributions, fullfile(outputDir, [tag '_distributions.png']), 'Resolution', 150);
    names = {'C', 'S', 'WC', 'WS', 'W'};
    stats = debug.stats;
    records = [stats.C; stats.S; stats.WC; stats.WS; stats.W];
    T = addvars(struct2table(records), string(names(:)), 'Before', 1, 'NewVariableNames', 'metric');
    T.validCount = repmat(stats.validCount, 5, 1);
    T.totalCount = repmat(stats.totalCount, 5, 1);
    ratioNames = fieldnames(stats.ratios);
    for k = 1:numel(ratioNames)
        name = ratioNames{k};
        T.(['ratio' name]) = repmat(stats.ratios.(name), 5, 1);
    end
    writetable(T, fullfile(outputDir, [tag '_stats.csv']));
    paramsUsed = p;
    prefix = displayOpts.prefix;
    save(fullfile(outputDir, [tag '_params.mat']), 'paramsUsed', 'prefix');
end
end
