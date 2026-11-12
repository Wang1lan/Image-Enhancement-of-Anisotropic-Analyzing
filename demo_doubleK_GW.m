clc;
clear;
close all;

% 方案 C：降噪只引导梯度权重，底图和带通细节仍取 Iwork。
% 使用：先运行 demo_gradientWeight，将标定门限复制到下方 GLow/GHigh。
projectDir = fileparts(mfilename('fullpath'));
addpath(fullfile(projectDir, 'functions'));
imgPath = fullfile(projectDir, 'testImg', 'WLI_1.jpg');
assert(isfile(imgPath), 'demo_doubleK_GW:ImageMissing', ...
    '找不到输入图像 %s，请检查 testImg 目录或修改 imgPath。', imgPath);
img = imread(imgPath);

enhanceMode = 'Y'; % 'RGB'、'Y' 或 'BOTH'；BOTH 为两条独立增强分支

% 主路径降噪：严格方案 C 保持 false；启用后必须在相同节点重新标定。
denoise.enable = false;
denoise.method = 'bilateral'; % 'bilateral' 或 'guided'
denoise.mode = 'Y'; % 'RGB' 联合三通道或仅 'Y' 降噪
denoise.windowSize = 5;
denoise.sigmaSpatial = 1.5;
denoise.sigmaRange = 10; % 原生 8 位灰度级
denoise.guidedEps = 4; % 原生灰度级平方

% 两频段、RGB 三通道共享同一亮度权重；门限不可沿用结构张量的 S 门限。
gradientWeight.enable = true; % false：NW 无自适应权重基线
gradientWeight.opts = struct('inputMax', 255, ...
    'GLow', 0.0024894783978878089, 'GHigh', 0.0071887207526713029, 'Wmin', 0.3, ...
    'validMask', [], 'debugLevel', 'off'); % 'off' / 'stats' / 'maps'
gradientWeight.opts.denoise = struct('enable', true, ...
    'method', 'bilateral', 'windowSize', 5, ...
    'sigmaSpatial', 1.5, 'sigmaRange', 10, 'guidedEps', 4);
gradientWeight.opts.mapping = struct('type', 'smoothstep', 'gamma', 2);
% type 为 'linear' / 'power' / 'smoothstep'；gamma 仅用于 power。
weightSmooth.windowSize = 9; % 映射后的权重均值滤波，边界为 symmetric
% 权重支路 denoise.enable=false：GW_raw 对照；true：方案 C。
% 输入节点、引导方法/参数改变后，应重新标定，不能静默沿用旧门限。
% 同一来源比较三种映射时，使用同一组 GLow/GHigh。

detailSign = 1;
% bandpassFilter_B 的 scaleLowFreq 固定为 1；这里只列出实际生效的核参数。
bp.fRadius_1 = 5;
bp.sigmaHighFreq_1 = 0.3;
bp.sigmaLowFreq_1 = 0.15;
detailGain_1 = 8;

bp.fRadius_2 = 3;
bp.sigmaHighFreq_2 = 0.55;
bp.sigmaLowFreq_2 = 0.3;
detailGain_2 = 4;

maskParams.Te = round(255 * 0.93);
maskParams.Th = round(255 * 0.85);
maskParams.Tl = round(255 * 0.10);
maskParams.dilateSizes = 5;
maskParams.numGaussian = 1;
evalParams.highlightDistance = 7;
evalParams.borderWidth = 2;

enableNoise = false;
noiseROI = []; % []：原图框选；也可填 [x y width height]，与标定 ROI 格式不同
noiseParams.patchSize = 5;
noiseParams.stride = 3;
noiseParams.tailNum = 5;
detailDisplayLimit = 30; % 权重和增益之前的带符号细节，单位为 8 位灰度级
saveResults = true;
outputDir = fullfile(projectDir, 'result', 'gradientWeight');

%% 必要初始化检查；关闭 GW 时不检查占位门限，也不运行引导分析
assert(isa(img, 'uint8') && ndims(img) == 3 && size(img,3) == 3, ...
    'demo_doubleK_GW:Input', '本增强 demo 的输入必须为 RGB uint8。');
assert(any(strcmp(enhanceMode, {'RGB','Y','BOTH'})), ...
    'demo_doubleK_GW:Mode', 'enhanceMode 必须为 RGB、Y 或 BOTH。');
if gradientWeight.enable
    assert(~isempty(gradientWeight.opts.GLow) && ~isempty(gradientWeight.opts.GHigh), ...
        'demo_doubleK_GW:ThresholdsRequired', ...
        '请先运行 demo_gradientWeight，将标定的 GLow/GHigh 填入本脚本顶部。');
    assert(gradientWeight.opts.inputMax == 255, 'demo_doubleK_GW:InputMax', ...
        '当前增强主流程固定为 8 位；其他满量程仅用于核心权重函数的独立验证。');
end
if denoise.enable
    assert(any(strcmp(denoise.method, {'bilateral','guided'})), ...
        'denoise.method 必须为 bilateral 或 guided。');
    assert(any(strcmp(denoise.mode, {'RGB','Y'})), 'denoise.mode 必须为 RGB 或 Y。');
    validateattributes(denoise.windowSize, {'numeric'}, ...
        {'real','scalar','finite','integer','positive','odd'});
    if strcmp(denoise.method, 'bilateral')
        validateattributes(denoise.sigmaSpatial, {'numeric'}, {'real','scalar','finite','positive'});
        validateattributes(denoise.sigmaRange, {'numeric'}, {'real','scalar','finite','positive'});
    else
        validateattributes(denoise.guidedEps, {'numeric'}, {'real','scalar','finite','positive'});
        assert(denoise.windowSize <= min(size(img,1), size(img,2)), ...
            '引导滤波窗口不能超过图像尺寸。');
        assert(exist('imguidedfilter', 'file') == 2, '主路径引导降噪需要 imguidedfilter。');
    end
end
runRGB = strcmp(enhanceMode, 'RGB') || strcmp(enhanceMode, 'BOTH');
runY = strcmp(enhanceMode, 'Y') || strcmp(enhanceMode, 'BOTH');
kernel_1 = bandpassFilter_B(bp.fRadius_1, bp.sigmaHighFreq_1, bp.sigmaLowFreq_1);
kernel_2 = bandpassFilter_B(bp.fRadius_2, bp.sigmaHighFreq_2, bp.sigmaLowFreq_2);

%% 原图亮斑掩模与共享增强输入；主路径降噪语义沿用 SW 版
[softMask, hardMask] = buildHighlightMask(img, maskParams);
I0 = double(img);
Iwork = I0;
denoiseLabel = '主路径降噪关闭';
analysisNode = 'I0 -> Y';
if denoise.enable
    if strcmp(denoise.method, 'bilateral')
        denoiseFcn = @(src) bilateralFilterShiftSum(src, denoise.windowSize, ...
            denoise.sigmaSpatial, denoise.sigmaRange);
    else
        denoiseFcn = @(src) imguidedfilter(src, src, ...
            'NeighborhoodSize', [denoise.windowSize denoise.windowSize], ...
            'DegreeOfSmoothing', denoise.guidedEps);
    end
    if strcmp(denoise.mode, 'RGB')
        Iwork = denoiseFcn(I0);
    else
        ycc0 = rgbToYCbCrFull(I0);
        yccDenoised = ycc0;
        yccDenoised(:,:,1) = denoiseFcn(ycc0(:,:,1));
        Iwork = yCbCrFullToRgb(yccDenoised);
    end
    denoiseLabel = sprintf('主路径 %s %s 降噪', denoise.mode, denoise.method);
    analysisNode = sprintf('Iwork -> Y (main %s %s)', denoise.method, denoise.mode);
end

%% 只在亮度分析支路降噪；Yguide 不回写 Iwork、Ywork 或两路细节
if runY || gradientWeight.enable
    yccWork = rgbToYCbCrFull(Iwork);
    Ywork = yccWork(:,:,1);
end
W = ones(size(softMask));
gradientDebug = struct();
weightLabel = 'NW';
source = 'off';
sourceLabel = 'off';
mappingLabel = 'off';
if gradientWeight.enable
    Yanalysis = min(max(Ywork, 0), 255); % 仅修正分析输入的浮点越界
    [W, gradientDebug] = computeGradientWeight(Yanalysis, gradientWeight.opts);
    source = 'raw';
    sourceLabel = 'raw';
    if gradientWeight.opts.denoise.enable && gradientWeight.opts.denoise.windowSize > 1
        source = 'denoised';
        sourceLabel = validatestring(gradientWeight.opts.denoise.method, {'bilateral','guided'});
    end
    mappingType = validatestring(gradientWeight.opts.mapping.type, {'linear','power','smoothstep'});
    mappingLabel = [upper(mappingType(1)) mappingType(2:end)];
    if strcmp(mappingType, 'power')
        mappingLabel = sprintf('Power_g%.17g', gradientWeight.opts.mapping.gamma);
    end
    weightLabel = ['GW_' sourceLabel '_' mappingLabel];
    fprintf('%s | source=%s | node=%s | mapping=%s, gamma=%.17g\n', ...
        weightLabel, source, analysisNode, mappingType, gradientWeight.opts.mapping.gamma);
    fprintf('GLow=%.17g, GHigh=%.17g, Wmin=%g, weightSmooth=%dx%d (symmetric)\n', ...
        gradientWeight.opts.GLow, gradientWeight.opts.GHigh, gradientWeight.opts.Wmin, ...
        weightSmooth.windowSize, weightSmooth.windowSize);
    fprintf('改变梯度来源、分析节点或权重支路降噪参数后，请重新标定门限。\n');
    if denoise.enable
        fprintf('组合实验：主路径降噪 + 梯度权重；标定必须使用相同 Iwork 节点和配置。\n');
    end
end
gradientWeight.sourceLabel = sourceLabel;
gradientWeight.mappingLabel = mappingLabel;

W_smooth = imboxfilt(W, weightSmooth.windowSize, 'Padding', 'symmetric');

detailWeight = (1-softMask) .* W_smooth;
if isfield(gradientDebug, 'W')
    gradientDebug.W_smooth = W_smooth;
    gradientDebug.detailWeight = detailWeight;
    gradientDebug.weightSmooth = weightSmooth;
    if ~saveResults
        showGradientWeightDebug(gradientDebug, struct('save', false));
    end
end
experimentLabel = [weightLabel ' | ' denoiseLabel];

%% 两条独立增强分支：底图与两路细节均来自 Iwork
if runRGB
    detailRGB_1 = detailSign .* filterShiftSum(Iwork, kernel_1, 'replicate');
    detailRGB_2 = detailSign .* filterShiftSum(Iwork, kernel_2, 'replicate');
    imgRawRGB = Iwork + detailGain_1 .* detailWeight .* detailRGB_1 ...
                     + detailGain_2 .* detailWeight .* detailRGB_2;
    imgEnhRGB = uint8(round(min(max(imgRawRGB, 0), 255)));
    tegRGB = evaluateTEGGlobal(img, imgEnhRGB, hardMask, evalParams);
    fprintf('[RGB | %s] TEGGlobal=%.9g\n', experimentLabel, tegRGB);
end
if runY
    detailY_1 = detailSign .* filterShiftSum(Ywork, kernel_1, 'replicate');
    detailY_2 = detailSign .* filterShiftSum(Ywork, kernel_2, 'replicate');
    yccEnh = yccWork; % 仅更新 Y，Cb/Cr 保持分支输入值
    yccEnh(:,:,1) = Ywork + detailGain_1 .* detailWeight .* detailY_1 ...
                         + detailGain_2 .* detailWeight .* detailY_2;
    imgRawY = yCbCrFullToRgb(yccEnh);
    imgEnhY = uint8(round(min(max(imgRawY, 0), 255)));
    tegY = evaluateTEGGlobal(img, imgEnhY, hardMask, evalParams);
    fprintf('[YCbCr-Y | %s] TEGGlobal=%.9g\n', experimentLabel, tegY);
end

%% 整图对比；保留原图 axes 供 PCA 框选
originalFigure = figure('Name', ['双带通增强对比 | ' experimentLabel], ...
    'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0.05 0.15 0.9 0.65]);
imageLayout = tiledlayout(originalFigure, 1, 1 + runRGB + runY, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
originalAxes = nexttile(imageLayout);
imshow(img, 'Parent', originalAxes);
title(originalAxes, '原图');
detailMaps = cell(0, 2);
branchNames = {};
if runRGB
    resultAxes = nexttile(imageLayout);
    imshow(imgEnhRGB, 'Parent', resultAxes);
    title(resultAxes, ['RGB 增强 | ' experimentLabel], 'Interpreter', 'none');
    Kr = 0.299; Kb = 0.114; Kg = 1 - Kr - Kb;
    detailMaps(end+1,:) = { ...
        Kr*detailRGB_1(:,:,1) + Kg*detailRGB_1(:,:,2) + Kb*detailRGB_1(:,:,3), ...
        Kr*detailRGB_2(:,:,1) + Kg*detailRGB_2(:,:,2) + Kb*detailRGB_2(:,:,3)};
    branchNames{end+1} = 'RGB';
end
if runY
    resultAxes = nexttile(imageLayout);
    imshow(imgEnhY, 'Parent', resultAxes);
    title(resultAxes, ['YCbCr-Y 增强 | ' experimentLabel], 'Interpreter', 'none');
    detailMaps(end+1,:) = {detailY_1, detailY_2};
    branchNames{end+1} = 'YCbCr-Y';
end

%% 掩模、梯度权重及增益之前的带符号细节，所有分支共用色阶
detailFigure = figure('Name', ['权重和增益之前的带符号细节 | ' experimentLabel], ...
    'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
detailLayout = tiledlayout(detailFigure, numel(branchNames), 2, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
colorRamp = linspace(0, 1, 129)';
colormap(detailFigure, [colorRamp colorRamp ones(129,1); ...
    ones(128,1) flipud(colorRamp(1:128)) flipud(colorRamp(1:128))]);
for branchIndex = 1:numel(branchNames)
    for kernelIndex = 1:2
        detailAxes = nexttile(detailLayout);
        imagesc(detailAxes, detailMaps{branchIndex, kernelIndex}, ...
            [-detailDisplayLimit detailDisplayLimit]);
        axis(detailAxes, 'image', 'off'); colorbar(detailAxes);
        title(detailAxes, sprintf('%s 核 %d | %s', ...
            branchNames{branchIndex}, kernelIndex, experimentLabel), 'Interpreter', 'none');
    end
end

%% 可选 PCA：仅在原图选坐标，不显示局部裁剪或放大图
if enableNoise
    if isempty(noiseROI)
        assert(exist('drawrectangle', 'file') == 2, '交互选区需要 drawrectangle。');
        figure(originalFigure);
        title(originalAxes, '框选平坦、无饱和亮斑区域，双击确认');
        roiHandle = drawrectangle(originalAxes);
        wait(roiHandle);
        assert(isvalid(roiHandle), '未确认噪声 ROI。');
        roiPosition = round(roiHandle.Position);
        roiLeft = max(1, roiPosition(1));
        roiTop = max(1, roiPosition(2));
        roiRight = min(size(img,2), roiPosition(1) + roiPosition(3) - 1);
        roiBottom = min(size(img,1), roiPosition(2) + roiPosition(4) - 1);
        noiseROI = [roiLeft roiTop roiRight-roiLeft+1 roiBottom-roiTop+1];
    end
    if runRGB
        noiseRGB = evaluateNoisePCA(img, imgEnhRGB, noiseROI, noiseParams);
        fprintf('[RGB] sigmaOrig=%.9g, sigmaEnh=%.9g, noiseGain=%.9g\n', ...
            noiseRGB.sigmaOrig, noiseRGB.sigmaEnh, noiseRGB.noiseGain);
    end
    if runY
        noiseY = evaluateNoisePCA(img, imgEnhY, noiseROI, noiseParams);
        fprintf('[YCbCr-Y] sigmaOrig=%.9g, sigmaEnh=%.9g, noiseGain=%.9g\n', ...
            noiseY.sigmaOrig, noiseY.sigmaEnh, noiseY.noiseGain);
    end
end

%% 按需保存；所有结果已经生成，参数只记录实际生效的三参数带通核
if saveResults
    if ~isfolder(outputDir)
        mkdir(outputDir);
    end
    [~, imgName] = fileparts(imgPath);
    prefix = sprintf('%s_%s_main%d-%s-%s', imgName, enhanceMode, denoise.enable, ...
        denoise.method, denoise.mode);
    tag = [prefix '_' weightLabel '_' char(datetime('now', 'Format', 'yyyyMMdd''T''HHmmssSSS'))];
    saveOption = '-v7';
    if isfield(gradientDebug, 'W')
        savedDebug = showGradientWeightDebug(gradientDebug, ...
            struct('save', true, 'outputDir', outputDir, 'prefix', prefix));
        tag = savedDebug.outputTag;
        saveOption = '-append'; % 将增强配置追加到含精确 G/W 的同一参数 MAT
    end
    if runRGB
        imwrite(imgEnhRGB, fullfile(outputDir, [tag '_RGB.png']));
    end
    if runY
        imwrite(imgEnhY, fullfile(outputDir, [tag '_Y.png']));
    end
    save(fullfile(outputDir, [tag '_params.mat']), 'imgPath', 'enhanceMode', ...
        'denoise', 'gradientWeight', 'analysisNode', 'bp', 'kernel_1', 'kernel_2', ...
        'detailSign', 'detailGain_1', 'detailGain_2', 'maskParams', 'evalParams', ...
        'source', 'sourceLabel', 'mappingLabel', 'weightSmooth', ...
        'W', 'W_smooth', 'detailWeight', saveOption);
    if isfield(gradientDebug, 'paramsUsed')
        paramsUsed = gradientDebug.paramsUsed;
        stats = gradientDebug.stats;
        save(fullfile(outputDir, [tag '_params.mat']), 'paramsUsed', 'stats', '-append');
    end
end
