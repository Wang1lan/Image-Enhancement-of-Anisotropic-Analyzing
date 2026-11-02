%% 结构权重标定：整图分析 -> 两类 ROI -> 候选门限 -> W
% 运行后，分别用鼠标左键点击每个矩形的两个对角点；Enter 可终止选取。
% 也可在下方填写 ROI=[x1 y1 x2 y2]（含端点），复现同一组像素。
% ROI 应避开视场边缘、反光和饱和区；默认不自动推断视场/高亮掩模。
% sigma/rho、输入节点或 Gamma 改变后，需要重新标定 S 门限。
% 仅使用基础 MATLAB；本例不生成工程带通细节或实际增强效果图。

clear W dbg roiStats;
scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, 'functions'));

%% 集中配置
imgPath = fullfile(scriptDir, 'testImg', 'WLI.jpg');
noiseROI = [];       % 平坦噪声区：[x1 y1 x2 y2]；空时鼠标选取
structureROI = [];   % 可信方向结构区，尽量减少平坦背景混入
saveResults = false;
outputDir = fullfile(scriptDir, 'result', 'structureWeight');
opts = struct('inputMax', 255, 'sigma', 1, 'rho', 2, 'CLow', .2, 'CHigh', .8, ...
    'SLow', 0, 'SHigh', 1, 'epsilon', 1e-12, 'validMask', [], 'debugLevel', 'maps');
% 上面的 SLow=0、SHigh=1 仅用于首次采集 C/S；占位 W 不参与增强判断。
% 更换输入文件时显式设置 inputMax；不要根据存储类型或当前最大值猜测。

%% 固定系数灰度分支；保留原始文件幅值，归一化由核心函数完成
img = imread(imgPath);
if size(img, 3) == 3
    Yanalysis = .299*double(img(:,:,1)) + .587*double(img(:,:,2)) + .114*double(img(:,:,3));
else
    Yanalysis = double(img);
end
analysisNode = 'raw-file gray (0.299R + 0.587G + 0.114B for RGB)';
[~, dbgCalib] = computeStructureWeight(Yanalysis, opts);
margin = ceil(3*opts.sigma) + 1 + ceil(3*opts.rho);
[~, imgName] = fileparts(imgPath);
fprintf('Input: %s | node: %s\n', imgName, analysisNode);
fprintf('inputMax=%g, sigma=%g, rho=%g; ROI boundary margin >= %d pixels.\n', ...
    opts.inputMax, opts.sigma, opts.rho, margin);

%% 在原图选 ROI；取样来自整图分析结果，不对裁剪图重新滤波
figROI = figure('Color', 'w', 'Name', '选择标定 ROI', 'NumberTitle', 'off');
axROI = axes('Parent', figROI);
if size(img, 3) == 3
    image(axROI, img);
else
    imagesc(axROI, img, [0 opts.inputMax]); colormap(axROI, gray(256));
end
axis(axROI, 'image'); axis(axROI, 'off'); hold(axROI, 'on');
[noiseROI, noiseMask] = chooseROI(axROI, noiseROI, '平坦噪声区', ...
    [0 .85 .85], dbgCalib.validMask, margin);
[structureROI, structureMask] = chooseROI(axROI, structureROI, '可信方向结构区', ...
    [1 .65 0], dbgCalib.validMask, margin);
title(axROI, sprintf('noise [%s] | structure [%s]', ...
    num2str(noiseROI), num2str(structureROI)), 'Interpreter', 'none');

noiseS = dbgCalib.S(noiseMask);
structureS = dbgCalib.S(structureMask);
candidateLow = samplePercentiles(noiseS, 95);
candidateHigh = samplePercentiles(structureS, 75);
fprintf('noise ROI [%s]: S-P95 = %.12g (SLow candidate)\n', num2str(noiseROI), candidateLow);
fprintf('structure ROI [%s]: S-P75 = %.12g (SHigh candidate)\n', num2str(structureROI), candidateHigh);

%% 即使候选门限重叠，也保留两类分布供重新选择 ROI
logNoise = log10(noiseS + opts.epsilon);
logStructure = log10(structureS + opts.epsilon);
limits = [min([logNoise; logStructure]), max([logNoise; logStructure])];
if limits(1) == limits(2)
    limits = limits + [-.5 .5];
end
edges = linspace(limits(1), limits(2), 65);
figROIS = figure('Color', 'w', 'Name', 'ROI 能量分布比较', 'NumberTitle', 'off');
axS = axes('Parent', figROIS);
histogram(axS, logNoise, edges, 'Normalization', 'probability', ...
    'FaceAlpha', .45, 'DisplayName', 'noise');
hold(axS, 'on');
histogram(axS, logStructure, edges, 'Normalization', 'probability', ...
    'FaceAlpha', .45, 'DisplayName', 'structure');
xline(axS, log10(candidateLow+opts.epsilon), '--r', 'noise P95', 'HandleVisibility', 'off');
xline(axS, log10(candidateHigh+opts.epsilon), '--k', 'structure P75', 'HandleVisibility', 'off');
xlabel(axS, 'log10(S + epsilon)'); ylabel(axS, 'ROI 内像素比例');
title(axS, 'ROI 能量分布（整图分析后取样）'); grid(axS, 'on'); legend(axS, 'Location', 'best');

if candidateHigh <= candidateLow
    error('demo_structureWeight:ThresholdOverlap', ...
        ['SHigh candidate <= SLow candidate: ROI distributions overlap. ' ...
        'Recheck ROIs, scales and image quality; no final W was generated.']);
end

%% 用真实候选门限重算；本次标定只适用于当前输入节点、尺度及标度
opts.SLow = candidateLow;
opts.SHigh = candidateHigh;
clear dbgCalib;
[W, dbg] = computeStructureWeight(Yanalysis, opts);
roiStats = [roiStatistics(dbg, noiseMask, 'noise'); ...
    roiStatistics(dbg, structureMask, 'structure')];
disp(roiStats);
prefix = sprintf('%s_noise_x%d-%d_y%d-%d_structure_x%d-%d_y%d-%d', imgName, ...
    noiseROI([1 3 2 4]), structureROI([1 3 2 4]));
displayOpts = struct('save', saveResults, 'outputDir', outputDir, 'prefix', prefix);
figs = showStructureWeightDebug(dbg, displayOpts);
if displayOpts.save
    tag = sprintf('%s_s%g_r%g_C%g-%g_S%g-%g_%s', prefix, opts.sigma, opts.rho, ...
        opts.CLow, opts.CHigh, opts.SLow, opts.SHigh, datestr(now, 'yyyymmddTHHMMSSFFF'));
    exportgraphics(figROIS, fullfile(outputDir, [tag '_roiS.png']), 'Resolution', 150);
    writetable(roiStats, fullfile(outputDir, [tag '_roiStats.csv']));
    fprintf('Saved debug outputs to: %s\n', outputDir);
end

%% 增强主程序最小接入示例（以下变量须来自实际工程，故不在 demo 中执行）
% opts.debugLevel = 'off';
% W = computeStructureWeight(Yanalysis, opts);
% enhanced = double(base) + gain .* W .* double(detail);
% base/detail 必须同标度；detail 为有符号细节；调用方负责限幅和颜色重建。
% RGB 共用亮度权重时：enhanced = double(base) + gain.*bsxfun(@times, W, double(detail));

function [roi, mask] = chooseROI(ax, roi, label, color, validMask, margin)
if isempty(roi)
    axes(ax);
    title(ax, sprintf('%s：左键点两个对角；避开反光/视场边缘至少 %d 像素', label, margin));
    [x, y, button] = ginput(2);
    if numel(x) ~= 2 || any(button ~= 1)
        error('demo_structureWeight:ROICancelled', 'Two left-click corners are required.');
    end
    roi = [min(x), min(y), max(x), max(y)];
end
if ~isnumeric(roi) || ~isreal(roi) || numel(roi) ~= 4 || any(~isfinite(roi(:)))
    error('demo_structureWeight:ROI', 'ROI must contain [x1 y1 x2 y2].');
end
roi = round(roi(:).');
[h, w] = size(validMask);
if roi(1) > roi(3) || roi(2) > roi(4) || any(roi(1:2) <= margin) ...
        || roi(3) > w-margin || roi(4) > h-margin
    error('demo_structureWeight:ROIBoundary', ...
        'ROI must keep at least %d pixels from the image boundary; coordinates were not clipped.', margin);
end
support = validMask(roi(2)-margin:roi(4)+margin, roi(1)-margin:roi(3)+margin);
if ~all(support(:))
    error('demo_structureWeight:ROIMaskBoundary', ...
        'ROI filtering support crosses invalidMask pixels; choose an interior ROI.');
end
if min(roi(3:4)-roi(1:2)+1) <= 2*margin
    warning('demo_structureWeight:SmallROI', ...
        'ROI is small relative to the analysis support; calibration may be unreliable. ROI was not enlarged.');
end
mask = false(h, w);
mask(roi(2):roi(4), roi(1):roi(3)) = true;
rectangle('Parent', ax, 'Position', [roi(1:2), roi(3:4)-roi(1:2)], ...
    'EdgeColor', color, 'LineWidth', 1.5);
text(ax, roi(1), roi(2), label, 'Color', color, 'VerticalAlignment', 'bottom');
drawnow;
end

function T = roiStatistics(d, mask, label)
names = {'C', 'S', 'WC', 'WS', 'W'};
values = zeros(5, 11);
for k = 1:5
    v = d.(names{k})(mask);
    values(k,:) = [numel(v), min(v), max(v), mean(v), std(v,1), ...
        samplePercentiles(v, [5 50 75 90 95 99])];
end
T = array2table(values, 'VariableNames', ...
    {'count', 'min', 'max', 'mean', 'std', 'P05', 'P50', 'P75', 'P90', 'P95', 'P99'});
T = addvars(T, repmat(string(label), 5, 1), string(names(:)), ...
    'Before', 1, 'NewVariableNames', {'region', 'metric'});
end

function q = samplePercentiles(v, percentages)
v = sort(v(:)).';
pos = 1+(numel(v)-1)*percentages/100;
q = v(floor(pos)) + (pos-floor(pos)).*(v(ceil(pos))-v(floor(pos)));
end
