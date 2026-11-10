%% 方案 C 梯度门限标定：整图分析 -> 两类 ROI -> 固定门限 -> W
% 使用：手填 ROI，或分别用鼠标左键点击矩形两角；Enter 可取消。
% ROI=[x1 y1 x2 y2]，一基且含端点；不同于 PCA 的 [x y width height]。
% 将输出的 GLow/GHigh 复制到 demo_doubleK_GW 顶部；改变分析节点或引导参数后重标定。
% 默认不保存，不自动推断黑色视场；ROI 支撑必须避开边界、显式无效区和亮斑保护区。
clc;
clear;
close all;
scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, 'functions'));

%% 集中配置：图像、主路径与引导降噪须与增强 demo 一致
imgPath = fullfile(scriptDir, 'testImg', 'WLI.jpg');
noiseROI = []; % 平坦噪声区：[x1 y1 x2 y2]；空时鼠标选取
structureROI = []; % 可信边缘/血管结构区，避免大量平坦背景
noisePercentile = 90;
structurePercentile = 75;
saveResults = false;
outputDir = fullfile(scriptDir, 'result', 'gradientWeightCalibration');

denoise.enable = false; % 严格 C：主路径关闭；组合实验须同步增强 demo 的配置
denoise.method = 'bilateral';
denoise.mode = 'Y';
denoise.windowSize = 5;
denoise.sigmaSpatial = 1.5;
denoise.sigmaRange = 10;
denoise.guidedEps = 4;

opts = struct('inputMax', 255, 'GLow', 0, 'GHigh', 1, 'Wmin', 0.3, ...
    'validMask', [], 'debugLevel', 'maps');
opts.denoise = struct('enable', true, 'method', 'bilateral', ...
    'windowSize', 5, 'sigmaSpatial', 1.5, 'sigmaRange', 10, 'guidedEps', 4);
% 首次的 GLow=0、GHigh=1 仅采集梯度分布；占位 W 不用于增强或最终权重展示。

maskParams.Te = round(255 * 0.93);
maskParams.Th = round(255 * 0.85);
maskParams.Tl = round(255 * 0.10);
maskParams.dilateSizes = 5;
maskParams.numGaussian = 1;

%% 复现增强端 Iwork -> Y 的分析节点；不提前裁剪 ROI 再滤波
img = imread(imgPath);
assert(isa(img, 'uint8') && ndims(img) == 3 && size(img,3) == 3, ...
    'demo_gradientWeight:Input', '与增强 demo 一致，输入必须为 RGB uint8。');
assert(opts.inputMax == 255, 'demo_gradientWeight:InputMax', '当前主流程的满量程固定为 255。');
validateattributes(noisePercentile, {'numeric'}, {'real','scalar','finite','>=',0,'<=',100});
validateattributes(structurePercentile, {'numeric'}, {'real','scalar','finite','>=',0,'<=',100});
I0 = double(img);
Iwork = I0;
mainRadius = 0;
analysisNode = 'I0 -> Y';
if denoise.enable
    assert(any(strcmp(denoise.method, {'bilateral','guided'})), ...
        'denoise.method 必须为 bilateral 或 guided。');
    assert(any(strcmp(denoise.mode, {'RGB','Y'})), 'denoise.mode 必须为 RGB 或 Y。');
    validateattributes(denoise.windowSize, {'numeric'}, ...
        {'real','scalar','finite','integer','positive','odd'});
    mainRadius = floor(denoise.windowSize/2);
    if strcmp(denoise.method, 'bilateral')
        validateattributes(denoise.sigmaSpatial, {'numeric'}, {'real','scalar','finite','positive'});
        validateattributes(denoise.sigmaRange, {'numeric'}, {'real','scalar','finite','positive'});
        denoiseFcn = @(src) bilateralFilterShiftSum(src, denoise.windowSize, ...
            denoise.sigmaSpatial, denoise.sigmaRange);
    else
        validateattributes(denoise.guidedEps, {'numeric'}, {'real','scalar','finite','positive'});
        assert(denoise.windowSize <= min(size(img,1), size(img,2)), ...
            '引导滤波窗口不能超过图像尺寸。');
        assert(exist('imguidedfilter', 'file') == 2, '主路径引导降噪需要 imguidedfilter。');
        denoiseFcn = @(src) imguidedfilter(src, src, ...
            'NeighborhoodSize', [denoise.windowSize denoise.windowSize], ...
            'DegreeOfSmoothing', denoise.guidedEps);
        mainRadius = 2*mainRadius;
    end
    if strcmp(denoise.mode, 'RGB')
        Iwork = denoiseFcn(I0);
    else
        ycc0 = rgbToYCbCrFull(I0);
        yccDenoised = ycc0;
        yccDenoised(:,:,1) = denoiseFcn(ycc0(:,:,1));
        Iwork = yCbCrFullToRgb(yccDenoised);
    end
    analysisNode = sprintf('Iwork -> Y (main %s %s)', denoise.method, denoise.mode);
end
yccWork = rgbToYCbCrFull(Iwork);
Yanalysis = min(max(yccWork(:,:,1), 0), 255); % 仅分析输入限幅，不回写底图
[~, dbgCalib] = computeGradientWeight(Yanalysis, opts);
d = dbgCalib.paramsUsed.denoise;
guideRadius = 0;
method = 'raw';
if d.enable && d.windowSize > 1
    guideRadius = floor(d.windowSize/2);
    method = d.method;
    if strcmp(d.method, 'guided')
        guideRadius = 2*guideRadius; % guided 的统计与系数平均扩大有效支撑
    end
end
margin = mainRadius + guideRadius + 1; % 末尾 1 为 Sobel 半径
[softMask, ~] = buildHighlightMask(img, maskParams);
calibrationMask = dbgCalib.validMask & (softMask == 0);
[~, imgName] = fileparts(imgPath);
fprintf('Input: %s | node: %s | GW_%s\n', imgName, analysisNode, method);
fprintf('inputMax=255, Sobel /8, L1; ROI 支撑边距至少 %d 像素。\n', margin);

%% 在原图选 ROI，样本取自已完成整图分析的 G
figROI = figure('Color', 'w', 'Name', '梯度标定：选择 ROI', 'NumberTitle', 'off');
axROI = axes('Parent', figROI);
imshow(img, 'Parent', axROI); hold(axROI, 'on');
[noiseROI, noiseMask] = chooseROI(axROI, noiseROI, '平坦噪声区', ...
    [0 .85 .85], calibrationMask, margin);
[structureROI, structureMask] = chooseROI(axROI, structureROI, '可信结构区', ...
    [1 .65 0], calibrationMask, margin);
title(axROI, sprintf('noise [%s] | structure [%s]', ...
    num2str(noiseROI), num2str(structureROI)), 'Interpreter', 'none');
noiseG = dbgCalib.G(noiseMask);
structureG = dbgCalib.G(structureMask);
candidateLow = samplePercentiles(noiseG, noisePercentile);
candidateHigh = samplePercentiles(structureG, structurePercentile);
fprintf('noise ROI [%s]: G-P%g = %.17g (GLow candidate)\n', ...
    num2str(noiseROI), noisePercentile, candidateLow);
fprintf('structure ROI [%s]: G-P%g = %.17g (GHigh candidate)\n', ...
    num2str(structureROI), structurePercentile, candidateHigh);
roiStats = [roiStatistics(noiseG, 'noise', 'G'); ...
    roiStatistics(structureG, 'structure', 'G')];

%% 共享线性色阶显示分布；计算、百分位数与映射始终使用线性 G
limits = [min([noiseG; structureG]), max([noiseG; structureG])];
if limits(1) == limits(2)
    limits = limits + [-1 1]*max(1e-6, abs(limits(1))*.01); % 只扩展显示范围
end
edges = linspace(limits(1), limits(2), 65);
figROIG = figure('Color', 'w', 'Name', '梯度标定：ROI 分布', 'NumberTitle', 'off');
axG = axes('Parent', figROIG);
histogram(axG, noiseG, edges, 'Normalization', 'probability', ...
    'FaceAlpha', .45, 'DisplayName', 'noise');
hold(axG, 'on');
histogram(axG, structureG, edges, 'Normalization', 'probability', ...
    'FaceAlpha', .45, 'DisplayName', 'structure');
xline(axG, candidateLow, '--r', sprintf('noise P%g', noisePercentile), 'HandleVisibility', 'off');
xline(axG, candidateHigh, '--k', sprintf('structure P%g', structurePercentile), 'HandleVisibility', 'off');
xlabel(axG, '归一化 Sobel L1 梯度 G'); ylabel(axG, 'ROI 内像素比例');
title(axG, '整图分析后的 ROI 分布及候选门限');
grid(axG, 'on'); legend(axG, 'Location', 'best');

% 即使门限重叠，也保留 ROI、候选值和分布；不生成伪成功权重。
calibration = struct('inputPath', imgPath, 'imageSize', size(img), ...
    'analysisNode', analysisNode, 'mainDenoise', denoise, ...
    'paramsUsed', dbgCalib.paramsUsed, 'maskParams', maskParams, ...
    'roiFormat', '[x1 y1 x2 y2], 1-based inclusive', ...
    'noiseROI', noiseROI, 'structureROI', structureROI, 'margin', margin, ...
    'noisePercentile', noisePercentile, 'structurePercentile', structurePercentile, ...
    'candidateLow', candidateLow, 'candidateHigh', candidateHigh, ...
    'status', 'threshold_overlap', 'MATLABVersion', version);
tag = sprintf('%s_GW_%s_G%.6g-%.6g_W%g_%s', imgName, method, ...
    candidateLow, candidateHigh, opts.Wmin, char(datetime('now', 'Format', 'yyyyMMdd''T''HHmmssSSS')));
if saveResults
    if ~isfolder(outputDir)
        mkdir(outputDir);
    end
    exportgraphics(figROI, fullfile(outputDir, [tag '_roi.png']), 'Resolution', 150);
    exportgraphics(figROIG, fullfile(outputDir, [tag '_roiG.png']), 'Resolution', 150);
    writetable(roiStats, fullfile(outputDir, [tag '_roiStats.csv']));
    save(fullfile(outputDir, [tag '_params.mat']), 'calibration', 'roiStats', 'noiseG', 'structureG');
end
if candidateHigh <= candidateLow
    error('demo_gradientWeight:ThresholdOverlap', ...
        ['GHigh candidate <= GLow candidate。请检查 ROI、引导强度和图像质量；' ...
        '未生成最终 W，也未调整候选门限。']);
end

%% 从原分析输入重算，避免把 Yguide 再次送入降噪
opts.GLow = candidateLow;
opts.GHigh = candidateHigh;
clear dbgCalib;
[W, dbg] = computeGradientWeight(Yanalysis, opts);
roiStats = [roiStats; roiStatistics(W(noiseMask), 'noise', 'W'); ...
    roiStatistics(W(structureMask), 'structure', 'W')];
calibration.paramsUsed = dbg.paramsUsed;
calibration.status = 'ok';
disp(roiStats);
fprintf('\n复制到 demo_doubleK_GW（须保持同一分析节点与引导配置）：\n');
fprintf('gradientWeight.opts.GLow = %.17g;\n', opts.GLow);
fprintf('gradientWeight.opts.GHigh = %.17g;\n', opts.GHigh);
figs = showGradientWeightDebug(dbg, struct('save', false));
if saveResults
    exportgraphics(figs.maps, fullfile(outputDir, [tag '_maps.png']), 'Resolution', 150);
    writetable(roiStats, fullfile(outputDir, [tag '_roiStats.csv']));
    save(fullfile(outputDir, [tag '_params.mat']), 'calibration', 'roiStats', 'W', 'dbg');
    fprintf('Saved calibration outputs to: %s\n', outputDir);
end

function [roi, mask] = chooseROI(ax, roi, label, color, validMask, margin)
if isempty(roi)
    axes(ax);
    title(ax, sprintf('%s：左键点击两个对角；支撑边距 %d 像素，Enter 取消', label, margin));
    [x, y, button] = ginput(2);
    if numel(x) ~= 2 || any(button ~= 1)
        error('demo_gradientWeight:ROICancelled', '需要两个鼠标左键角点；标定已取消。');
    end
    roi = [min(x), min(y), max(x), max(y)];
end
validateattributes(roi, {'numeric'}, {'real','finite','vector','numel',4});
roi = round(roi(:).');
[h, w] = size(validMask);
assert(roi(1) <= roi(3) && roi(2) <= roi(4) && all(roi(1:2) > margin) && ...
    roi(3) <= w-margin && roi(4) <= h-margin, ...
    'demo_gradientWeight:ROIBoundary', ...
    'ROI 必须在图内并距边界至少 %d 像素；坐标未被平移或裁剪。', margin);
support = validMask(roi(2)-margin:roi(4)+margin, roi(1)-margin:roi(3)+margin);
assert(all(support(:)), 'demo_gradientWeight:ROIExcluded', ...
    'ROI 的滤波支撑跨越无效区或亮斑保护区，请重新选取。');
if min(roi(3:4)-roi(1:2)+1) <= 2*margin
    warning('demo_gradientWeight:SmallROI', 'ROI 相对滤波支撑较小；未自动扩大选区。');
end
mask = false(h, w);
mask(roi(2):roi(4), roi(1):roi(3)) = true;
rectangle('Parent', ax, 'Position', [roi(1:2), roi(3:4)-roi(1:2)], ...
    'EdgeColor', color, 'LineWidth', 1.5);
text(ax, roi(1), roi(2), label, 'Color', color, 'VerticalAlignment', 'bottom');
drawnow;
end

function T = roiStatistics(v, region, metric)
values = [numel(v), min(v), max(v), mean(v), std(v,1), ...
    samplePercentiles(v, [5 50 75 90 95 99])];
T = array2table(values, 'VariableNames', ...
    {'count','min','max','mean','std','P05','P50','P75','P90','P95','P99'});
T = addvars(T, string(region), string(metric), 'Before', 1, ...
    'NewVariableNames', {'region','metric'});
end

function q = samplePercentiles(v, percentages)
v = sort(v(:)).';
pos = 1+(numel(v)-1)*percentages/100;
q = v(floor(pos)) + (pos-floor(pos)).*(v(ceil(pos))-v(floor(pos)));
end
