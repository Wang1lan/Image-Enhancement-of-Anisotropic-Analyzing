clc;
clear;
close all;

% 双带通增强实验。参数集中在下方，输入为 RGB uint8。
projectDir = fileparts(mfilename('fullpath'));
imgPath = fullfile(projectDir, 'testImg', 'WLI.jpg');

enhanceMode = 'Y'; % 'RGB'、'Y' 或 'BOTH'；BOTH 独立对比两种增强
denoise.enable = false; % 与 PCA 评价开关 enableNoise 独立
denoise.method = 'bilateral'; % 'bilateral' 双边或 'guided' 引导去噪
denoise.mode = 'Y'; % 'RGB' 三通道或仅 'Y' 去噪，与增强模式独立
denoise.windowSize = 5; % 实际窗口边长，正奇数；1 表示不做邻域平均
denoise.sigmaSpatial = 1.5; % 仅双边，单位：像素
denoise.sigmaRange = 10; % 仅双边，单位：8 位灰度级
denoise.guidedEps = 4; % 仅引导，单位：灰度级平方

detailSign = 1; % 默认核主通带已用带符号响应和正弦输入核对为同相

% 低中频段
bp.fRadius_1 = 5;
bp.sigmaHighFreq_1 = 0.3;
bp.sigmaLowFreq_1 = 0.15;
bp.scaleLowFreq_1 = 0.8;      
detailGain_1 = 3;     % 叠加到原图的增量系数

% 中高频段
bp.fRadius_2 = 3;
bp.sigmaHighFreq_2 = 55;
bp.sigmaLowFreq_2 = 0.3;
bp.scaleLowFreq_2 = 0.95;
detailGain_2 = 2;

% 两路细节共用 Iwork 亮度生成的结构权重；先标定 S 门限，再启用。
structureWeight.enable = false;
structureWeight.opts = struct('inputMax', 255, 'sigma', 1, 'rho', 2, ...
    'CLow', .2, 'CHigh', .8, 'SLow', [], 'SHigh', [], 'epsilon', 1e-12, ...
    'validMask', [], 'debugLevel', 'off'); % 'maps' 显示调试图，'stats' 仅保留统计
% S 门限应在同一 Iwork 输入节点、sigma/rho 下标定；不要沿用其他流程的门限。

% 手动启用保存示例；分支、去噪方法/通道和核参数区分结果，使用无损 PNG。
%{
savePath = fullfile(projectDir, 'result');
saveName = ['-', num2str(bp.fRadius_1), '_', ...
            num2str(bp.sigmaHighFreq_1), '_', num2str(bp.sigmaLowFreq_1),...
            '_', num2str(bp.scaleLowFreq_1), '_', num2str(detailGain_1), '-', num2str(bp.fRadius_2), '_', ...
            num2str(bp.sigmaHighFreq_2), '_', num2str(bp.sigmaLowFreq_2),...
            '_', num2str(bp.scaleLowFreq_2), '_', num2str(detailGain_2)];
if denoise.enable == true
    isdeniose = ['-dN-' denoise.method '-' denoise.mode];
elseif denoise.enable == false
    isdeniose = '-N';
end
if runRGB
    imwrite(imgEnhRGB, fullfile(savePath, ['RGB-' isdeniose saveName '.png']));
end
if runY
    imwrite(imgEnhY, fullfile(savePath, ['Y-' isdeniose saveName '.png']));
end
%}

maskParams.Te = round(255 * 0.93);
maskParams.Th = round(255 * 0.85);
maskParams.Tl = round(255 * 0.10);
maskParams.dilateSizes = 5;
maskParams.numGaussian = 1;

evalParams.highlightDistance = 7;
evalParams.borderWidth = 2;

enableNoise = false;
noiseROI = [];      % []：在原图框选并双击确认；也可填 [x y width height]
noiseParams.patchSize = 5;
noiseParams.stride = 3;
noiseParams.tailNum = 5;

detailDisplayLimit = 30; % 带符号细节热图范围 [-L,L]，单位为 8 位灰度级

%% 资源与必要初始化检查
assert(isfolder(fullfile(projectDir, 'Functions')) && isfile(imgPath), ...
    '缺少 Functions 目录或输入图像。');
addpath(fullfile(projectDir, 'Functions'));
assert(exist('fsamp2', 'file') == 2 && exist('padarray', 'file') == 2 && ...
    exist('bwdist', 'file') == 2, '需要 Image Processing Toolbox。');
img = imread(imgPath);
assert(isa(img, 'uint8') && ndims(img) == 3 && size(img, 3) == 3, ...
    '输入应为 RGB uint8 图像。');
assert(isscalar(bp.fRadius_1) && isfinite(bp.fRadius_1) && bp.fRadius_1 >= 1 && ...
    bp.fRadius_1 == floor(bp.fRadius_1) && ...
    all(isfinite([bp.sigmaHighFreq_1 bp.sigmaLowFreq_1 bp.scaleLowFreq_1])) && ...
    bp.sigmaHighFreq_1 > bp.sigmaLowFreq_1 && bp.sigmaLowFreq_1 > 0 && ...
    bp.scaleLowFreq_1 >= 0, '请检查核半径和频域响应参数。');
assert(isscalar(bp.fRadius_2) && isfinite(bp.fRadius_2) && bp.fRadius_2 >= 1 && ...
    bp.fRadius_2 == floor(bp.fRadius_2) && ...
    all(isfinite([bp.sigmaHighFreq_2 bp.sigmaLowFreq_2 bp.scaleLowFreq_2])) && ...
    bp.sigmaHighFreq_2 > bp.sigmaLowFreq_2 && bp.sigmaLowFreq_2 > 0 && ...
    bp.scaleLowFreq_2 >= 0, '请检查核半径和频域响应参数。');
assert(isscalar(detailSign) && ismember(detailSign, [-1 1]) && ...
    isscalar(detailGain_1) && isfinite(detailGain_1) && ...
    isscalar(detailGain_2) && isfinite(detailGain_2) && ...
    isscalar(detailDisplayLimit) && isfinite(detailDisplayLimit) && ...
    detailDisplayLimit > 0, '请检查细节符号、增量增益和热图范围。');
assert(all(isfinite([maskParams.Te maskParams.Th maskParams.Tl])) && ...
    0 <= maskParams.Th && maskParams.Th <= maskParams.Te && ...
    maskParams.Te <= 255 && 0 <= maskParams.Tl && maskParams.Tl <= 255 && ...
    all(isfinite(maskParams.dilateSizes)) && ...
    all(maskParams.dilateSizes >= 1 & mod(maskParams.dilateSizes, 2) == 1) && ...
    isscalar(maskParams.numGaussian) && isfinite(maskParams.numGaussian) && ...
    maskParams.numGaussian >= 0 && ...
    maskParams.numGaussian == floor(maskParams.numGaussian), '请检查掩模参数。');
assert(all(isfinite([evalParams.highlightDistance evalParams.borderWidth])) && ...
    evalParams.highlightDistance >= 0 && evalParams.borderWidth >= 0 && ...
    evalParams.borderWidth == floor(evalParams.borderWidth), '请检查评价区域参数。');
assert(any(strcmp(enhanceMode, {'RGB', 'Y', 'BOTH'})), ...
    'enhanceMode 必须为 RGB、Y 或 BOTH。');
if denoise.enable
    assert(any(strcmp(denoise.method, {'bilateral', 'guided'})), ...
        'denoise.method 必须为 bilateral 或 guided。');
    assert(any(strcmp(denoise.mode, {'RGB', 'Y'})), 'denoise.mode 必须为 RGB 或 Y。');
    assert(isscalar(denoise.windowSize) && isfinite(denoise.windowSize) && ...
        denoise.windowSize > 0 && mod(denoise.windowSize, 2) == 1, ...
        '去噪窗口边长必须为有限正奇整数。');
    switch denoise.method
        case 'bilateral'
            assert(isscalar(denoise.sigmaSpatial) && isfinite(denoise.sigmaSpatial) && ...
                denoise.sigmaSpatial > 0 && isscalar(denoise.sigmaRange) && ...
                isfinite(denoise.sigmaRange) && denoise.sigmaRange > 0, ...
                '去噪 sigma 必须为有限正数。');
        case 'guided'
            assert(isscalar(denoise.guidedEps) && isfinite(denoise.guidedEps) && ...
                denoise.guidedEps > 0, 'guidedEps 必须为有限正标量。');
            assert(denoise.windowSize <= min(size(img,1), size(img,2)), ...
                '引导滤波窗口边长不能超过图像尺寸。');
            assert(exist('imguidedfilter', 'file') == 2, '引导去噪需要 imguidedfilter。');
    end
end
runRGB = strcmp(enhanceMode, 'RGB') || strcmp(enhanceMode, 'BOTH');
runY = strcmp(enhanceMode, 'Y') || strcmp(enhanceMode, 'BOTH');

kernel_1 = bandpassFilter_B(bp.fRadius_1, bp.sigmaHighFreq_1, ...
    bp.sigmaLowFreq_1);
assert(isequal(size(kernel_1), [2*bp.fRadius_1+1 2*bp.fRadius_1+1]) && ...
    all(isfinite(kernel_1(:))) && any(kernel_1(:) > 0) && any(kernel_1(:) < 0) && ...
    abs(sum(kernel_1(:))) <= 1e-10 * sum(abs(kernel_1(:))), ...
    '带通核尺寸、有限性或零和检查失败，请检查归一化参数。');

kernel_2 = bandpassFilter_B(bp.fRadius_2, bp.sigmaHighFreq_2, ...
    bp.sigmaLowFreq_2);
assert(isequal(size(kernel_2), [2*bp.fRadius_2+1 2*bp.fRadius_2+1]) && ...
    all(isfinite(kernel_2(:))) && any(kernel_2(:) > 0) && any(kernel_2(:) < 0) && ...
    abs(sum(kernel_2(:))) <= 1e-10 * sum(abs(kernel_2(:))), ...
    '带通核尺寸、有限性或零和检查失败，请检查归一化参数。');

%% 原图掩模与共享去噪输入
[softMask, hardMask] = buildHighlightMask(img, maskParams);
I0 = double(img);
Iwork = I0;
denoiseLabel = '去噪关闭';
if denoise.enable
    switch denoise.method
        case 'bilateral'
            denoiseFcn = @(src) bilateralFilterShiftSum(src, denoise.windowSize, ...
                denoise.sigmaSpatial, denoise.sigmaRange);
            denoiseLabel = [denoise.mode ' 双边去噪'];
            fprintf('method=%s, mode=%s, windowSize=%d, sigmaSpatial=%.9g, sigmaRange=%.9g\n', ...
                denoise.method, denoise.mode, denoise.windowSize, ...
                denoise.sigmaSpatial, denoise.sigmaRange);
        case 'guided'
            denoiseFcn = @(src) imguidedfilter(src, src, ...
                'NeighborhoodSize', [denoise.windowSize denoise.windowSize], ...
                'DegreeOfSmoothing', denoise.guidedEps);
            denoiseLabel = [denoise.mode ' 引导去噪'];
            fprintf('method=%s, mode=%s, windowSize=%d, guidedEps=%.9g\n', ...
                denoise.method, denoise.mode, denoise.windowSize, denoise.guidedEps);
    end
    if strcmp(denoise.mode, 'RGB')
        Iwork = denoiseFcn(I0);
    else
        ycc0 = rgbToYCbCrFull(I0);
        yccDenoised = ycc0;
        yccDenoised(:,:,1) = denoiseFcn(ycc0(:,:,1));
        Iwork = yCbCrFullToRgb(yccDenoised);
    end
end

%% 共享亮度分析与融合权重；结构分析不改变细节提取输入
if runY || structureWeight.enable
    yccWork = rgbToYCbCrFull(Iwork);
    Ywork = yccWork(:,:,1);
end
detailWeight = 1 - softMask;
if structureWeight.enable
    % 8 位亮度理论范围为 [0,255]，仅在分析分支修正去噪/颜色往返的浮点越界。
    Yanalysis = min(max(Ywork, 0), 255);
    [W, structureDebug] = computeStructureWeight(Yanalysis, structureWeight.opts);
    detailWeight = detailWeight .* W;
    if strcmp(structureWeight.opts.debugLevel, 'maps')
        showStructureWeightDebug(structureDebug, struct('save', false));
    end
end

%% 两条独立增强分支：底图与两路细节均来自 Iwork
if runRGB
    detailRGB_1 = detailSign .* filterShiftSum(Iwork, kernel_1, 'replicate');
    detailRGB_2 = detailSign .* filterShiftSum(Iwork, kernel_2, 'replicate');
    imgRawRGB = Iwork + detailGain_1 .* detailWeight .* detailRGB_1 ...
                     + detailGain_2 .* detailWeight .* detailRGB_2;
    imgEnhRGB = uint8(round(min(max(imgRawRGB, 0), 255)));
    tegRGB = evaluateTEGGlobal(img, imgEnhRGB, hardMask, evalParams);
    fprintf('[RGB] TEGGlobal=%.9g\n', tegRGB);
end
if runY
    detailY_1 = detailSign .* filterShiftSum(Ywork, kernel_1, 'replicate');
    detailY_2 = detailSign .* filterShiftSum(Ywork, kernel_2, 'replicate');
    yccEnh = yccWork; % 仅更新 Y，Cb、Cr 保持分支输入值
    yccEnh(:,:,1) = Ywork + detailGain_1 .* detailWeight .* detailY_1 ...
                         + detailGain_2 .* detailWeight .* detailY_2;
    imgRawY = yCbCrFullToRgb(yccEnh);
    imgEnhY = uint8(round(min(max(imgRawY, 0), 255)));
    tegY = evaluateTEGGlobal(img, imgEnhY, hardMask, evalParams);
    fprintf('[YCbCr-Y] TEGGlobal=%.9g\n', tegY);
end

%% 整图对比；保留原图 axes 供 PCA 框选
originalFigure = figure('Name', ['双带通增强对比 | ' denoiseLabel], ...
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
    title(resultAxes, ['RGB 增强 | ' denoiseLabel]);
    Kr = 0.299; Kb = 0.114; Kg = 1 - Kr - Kb;
    detailMaps(end+1,:) = { ...
        Kr*detailRGB_1(:,:,1) + Kg*detailRGB_1(:,:,2) + Kb*detailRGB_1(:,:,3), ...
        Kr*detailRGB_2(:,:,1) + Kg*detailRGB_2(:,:,2) + Kb*detailRGB_2(:,:,3)};
    branchNames{end+1} = 'RGB';
end
if runY
    resultAxes = nexttile(imageLayout);
    imshow(imgEnhY, 'Parent', resultAxes);
    title(resultAxes, ['YCbCr-Y 增强 | ' denoiseLabel]);
    detailMaps(end+1,:) = {detailY_1, detailY_2};
    branchNames{end+1} = 'YCbCr-Y';
end

%% 掩模、结构权重及增益之前的带符号细节，所有分支共用色阶
detailFigure = figure('Name', ['带符号细节 | ' denoiseLabel], ...
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
        axis(detailAxes, 'image', 'off');
        colorbar(detailAxes);
        title(detailAxes, sprintf('%s 核 %d | %s', ...
            branchNames{branchIndex}, kernelIndex, denoiseLabel));
    end
end
title(detailLayout, '掩模、结构权重和增益之前的带符号亮度细节');

% olymEnhance = imread(fullfile(projectDir, 'testImg', 'olymEnhance.jpg'));
% figure('Name', '奥巴增强图像', 'NumberTitle', 'off');
% imshow(olymEnhance);
% title('奥巴增强图像');

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
        % 取框选范围与图像的实际交集，不将越界矩形平移到图内。
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
