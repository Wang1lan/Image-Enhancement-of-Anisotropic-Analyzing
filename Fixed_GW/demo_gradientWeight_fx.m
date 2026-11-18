function report = demo_gradientWeight_fx()
%DEMO_GRADIENTWEIGHT_FX 验证12位梯度权重，显示并保存逐阶段比较。
%   在 Fixed_GW 目录运行 report = demo_gradientWeight_fx。
%   JPG 仅用于构造测试 Y12，所有验收矩阵保存在 MAT 文件中。

% 数值范围规定
% Wmin = [0.1, 0.5];
% GLow = [0, 0.1];
% GHigh = [0, 0.1];
% min(GHigh-GLow) = 0.003， 两阈值的最小间隔为0.003


fixedDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(fixedDir);
oldPath = path;
pathCleanup = onCleanup(@() path(oldPath));
addpath(fullfile(projectDir, 'functions'), fullfile(fixedDir, 'functions'));

GLow = 0.0024894783978878089;
GHigh = 0.0071887207526713029;
Wmin = 0.3;
saveResults = true;
outputDir = fullfile(fixedDir, 'result', 'gradientWeightFixed', 'Q6_Q17_7x7_COMEN_approx49');
imgPath = fullfile(projectDir, 'testImg', 'WLI_1.jpg');
prefix = 'WLI_1';
if saveResults && ~isfolder(outputDir)
    mkdir(outputDir);
end

testResults = runtests(fullfile(fixedDir, 'tests', 'test_gradientWeight_fx.m'));
testSummary = table(string({testResults.Name}).', [testResults.Passed].', ...
    [testResults.Failed].', [testResults.Incomplete].', [testResults.Duration].', ...
    'VariableNames', {'Name','Passed','Failed','Incomplete','DurationSeconds'});
if saveResults
    writetable(testSummary, fullfile(outputDir, [prefix '_tests.csv']));
end
assertSuccess(testResults);

rgb8 = imread(imgPath);
ycc01 = comen_rgb2ycbcr_ft(double(rgb8)/255);
Y12 = uint16(round(4095*ycc01(:,:,1)));
params = prepareGradientWeightFxParams(GLow, GHigh, Wmin);
regs = loadGradientWeightFxRegs(params.Wmin12,params.Tlow_Q6,params.Thigh_Q6, ...
    params.K_Q17);
[W12, debugFx] = computeGradientWeight_fx(Y12, regs);
[Wfloat, debugFloat] = computeGradientWeight_floatRef(Y12, GLow, GHigh, Wmin);

errors = {double(debugFx.G)-round(32760*debugFloat.G), ...
    double(debugFx.W0_12)-round(4095*debugFloat.W0), ...
    double(W12)-round(4095*Wfloat)};
boundary = false(size(Y12));
boundary(1:4,:) = true;
boundary(end-3:end,:) = true;
boundary(:,1:4) = true;
boundary(:,end-3:end) = true;
threshold = abs(double(debugFx.G_Q6)-double(params.Tlow_Q6)) <= 64 | ...
    abs(double(debugFx.G_Q6)-double(params.Thigh_Q6)) <= 64;
regions = {true(size(Y12)), boundary, ~boundary, threshold};
regionNames = ["All", "Boundary4", "Interior", "Threshold1Gradient"];
stageNames = ["G", "W0", "Wfinal"];
scales = [32760 4095 4095];
metrics = table;
for stage = 1:numel(stageNames)
    for region = 1:numel(regionNames)
        row = metricRow(stageNames(stage), regionNames(region), ...
            errors{stage}(regions{region}), scales(stage));
        metrics = [metrics; row]; %#ok<AGROW>
    end
end

mapping.G = uint16(0:24570);
mapping.W0Fixed = mapGradientWeight_fx(mapping.G, params);
mapping.W0Float = Wmin + (1-Wmin)*min(max( ...
    (double(mapping.G)/32760-GLow)/(GHigh-GLow), 0), 1);
mapping.ErrorLSB = double(mapping.W0Fixed)-round(4095*mapping.W0Float);
sobelPass = isequal(debugFx.Gx, int32(round(32760*debugFloat.Gx))) && ...
    isequal(debugFx.Gy, int32(round(32760*debugFloat.Gy))) && ...
    all(errors{1}(:) == 0);
rangePass = isa(W12, 'uint16') && isequal(size(W12), size(Y12)) && ...
    all(W12(:) >= params.Wmin12 & W12(:) <= uint16(4095));
report.unitTestsPass = all([testResults.Passed]);
report.integerValidationPass = report.unitTestsPass && sobelPass && rangePass;
report.mappingMaxErrorLSB = max(abs(mapping.ErrorLSB));
report.imageW0MaxErrorLSB = max(abs(errors{2}(:)));
report.imageFinalMaxErrorLSB = max(abs(errors{3}(:)));
report.accuracyPass = report.mappingMaxErrorLSB <= 1 && ...
    report.imageW0MaxErrorLSB <= 1 && report.imageFinalMaxErrorLSB <= 2;
report.metrics = metrics;
report.params = params;
report.colorSpace = 'COMEN RGB/255 to YCbCr, zero-centered chroma; Y12=round(4095*Y01)';
report.outputDir = outputDir;

fprintf('\nSource: %s\n', imgPath);
fprintf('Size: %d x %d, Y12 range: [%d,%d]\n', size(Y12,1), size(Y12,2), ...
    min(Y12(:)), max(Y12(:)));
fprintf('GLow=%.17g, GHigh=%.17g, Wmin=%.17g\n', GLow, GHigh, Wmin);
fprintf('Q6/%s  R0..R3=[%u %u %u %u]  D=%u  C=%u  s=%d\n', ...
    params.Format,regs.Wmin12,regs.Tlow_Q6,regs.Thigh_Q6,regs.SlopeWord,regs.D_Q6,regs.C,regs.s);
disp(metrics(:, {'Stage','Region','Count','MAE_LSB','RMSE_LSB','Max_LSB', ...
    'P95Abs_LSB','P99Abs_LSB','ExactRatio'}));
fprintf('Integer validation: %s, float accuracy targets: %s\n', ...
    passLabel(report.integerValidationPass), passLabel(report.accuracyPass));
fprintf('Mapping / image W0 / image Wfinal max error: %g / %g / %g LSB\n', ...
    report.mappingMaxErrorLSB, report.imageW0MaxErrorLSB, report.imageFinalMaxErrorLSB);

outliers = table;
if ~report.accuracyPass || ~sobelPass
    outliers = gradientWeightAccuracyOutliers(debugFx, debugFloat, regs);
    disp(outliers(1:min(20,height(outliers)),:));
    warning('demo_gradientWeight_fx:Accuracy', ...
        '浮点精度目标未全部达到，误差坐标及分解见 outliers 表与 MAT。');
end

figG = compareFigure(32760*debugFloat.G, double(debugFx.G), errors{1}, ...
    [0 max(1,max(double(debugFx.G(:))))], 'Sobel L1（原始梯度单位）');
figW0 = compareFigure(debugFloat.W0, double(debugFx.W0_12)/4095, errors{2}, ...
    [0 1], '初始权重');
figFinal = compareFigure(Wfloat, double(W12)/4095, errors{3}, [0 1], '7×7最终权重');
figMapping = figure('Name', 'Q6/Q17 mapping', 'Color', 'w', ...
    'Position', [120 120 1100 650]);
layout = tiledlayout(figMapping, 2, 1, 'TileSpacing', 'compact');
nexttile(layout);
plot(double(mapping.G), mapping.W0Float, 'LineWidth', 1.4);
hold on;
plot(double(mapping.G), double(mapping.W0Fixed)/4095, '--', 'LineWidth', 1.2);
xlim([0 ceil(double(params.Thigh_Q6)/64*1.3)]);
ylim([0 1]);
ylabel('权重');
legend('浮点线性映射', '定点权重 /4095', 'Location', 'southeast');
grid on;
nexttile(layout);
stairs(double(mapping.G), mapping.ErrorLSB, 'LineWidth', 1.1);
xlim([0 ceil(double(params.Thigh_Q6)/64*1.3)]);
ylim([-1.5 1.5]);
xlabel('Sobel L1 原始整数梯度');
ylabel('定点 − 浮点量化值 / LSB');
grid on;

if saveResults
    imwrite(double(Y12)/4095, fullfile(outputDir, [prefix '_Y12.png']));
    figures = [figG figW0 figFinal figMapping];
    suffixes = {'G_compare','W0_compare','Wfinal_compare','mapping_compare'};
    for k = 1:numel(figures)
        exportgraphics(figures(k), fullfile(outputDir, ...
            [prefix '_' suffixes{k} '.png']), 'Resolution', 140);
    end
    writetable(metrics, fullfile(outputDir, [prefix '_metrics.csv']));
    writetable(parameterTable(params), fullfile(outputDir, [prefix '_params.csv']));
    if ~isempty(outliers)
        writetable(outliers, fullfile(outputDir, [prefix '_outliers.csv']));
    end
    matlabVersion = version;
    save(fullfile(outputDir, [prefix '_debug.mat']), 'Y12', 'params', 'regs', ...
        'debugFx', 'debugFloat', 'metrics', 'mapping', 'testSummary', ...
        'report', 'outliers', 'imgPath', 'matlabVersion', '-v7');
    fprintf('Saved results: %s\n', outputDir);
end
assert(report.integerValidationPass, 'demo_gradientWeight_fx:IntegerValidation', ...
    '整数一致性、Sobel或输出范围验收失败。');
end

function row = metricRow(stage, region, errorValues, scale)
names = {'Count','MAE_LSB','RMSE_LSB','Max_LSB','P95Abs_LSB','P99Abs_LSB', ...
    'ExactRatio','Over1Ratio','Over2Ratio','Over4Ratio','Over8Ratio', ...
    'MAE_Normalized','RMSE_Normalized','Max_Normalized', ...
    'P95Abs_Normalized','P99Abs_Normalized'};
count = numel(errorValues);
values = [count nan(1,numel(names)-1)];
if count > 0
    absolute = sort(abs(errorValues(:)));
    positions = 1+(count-1)*[0.95 0.99];
    percentiles = absolute(floor(positions)).' + ...
        (positions-floor(positions)).*(absolute(ceil(positions)).'-absolute(floor(positions)).');
    magnitudes = [mean(absolute), sqrt(mean(errorValues(:).^2)), ...
        absolute(end), percentiles];
    ratios = [mean(absolute == 0), mean(absolute > 1), mean(absolute > 2), ...
        mean(absolute > 4), mean(absolute > 8)];
    values = [count magnitudes ratios magnitudes/scale];
end
row = [table(stage, region, 'VariableNames', {'Stage','Region'}), ...
    array2table(values, 'VariableNames', names)];
end

function fig = compareFigure(reference, fixed, errorValues, limits, caption)
fig = figure('Name', caption, 'Color', 'w', 'Position', [80 160 1500 480]);
layout = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile(layout);
imagesc(reference, limits);
axis image off;
title([caption ' · 浮点参考']);
colorbar;
nexttile(layout);
imagesc(fixed, limits);
axis image off;
title([caption ' · 整数模型']);
colorbar;
nexttile(layout);
imagesc(abs(errorValues), [0 max(1,max(abs(errorValues(:))))]);
axis image off;
title('相对浮点量化值的绝对误差 / LSB');
colorbar;
colormap(fig, parula);
end

function values = parameterTable(p)
Name = ["GLow";"GHigh";"Wmin";"Wmin12";"Tlow_Q6";"Thigh_Q6"; ...
    "SlopeWord";"D_Q6";"C";"s";"fracBitsGradient";"inputMax"];
Value = compose("%.17g",[p.GLow;p.GHigh;p.Wmin;double(p.Wmin12); ...
    double(p.Tlow_Q6);double(p.Thigh_Q6);double(p.SlopeWord);double(p.D_Q6); ...
    double(p.C);p.s;6;4095]);
DeviceFormat = repmat(p.Format,size(Name));
Format = ["double";"double";"double";"u12";"u18 Q6";"u18 Q6"; ...
    "u17 Q17";"u18 Q6";"u17";"control";"control";"u12"];
Unit = ["normalized Sobel L1";"normalized Sobel L1";"weight"; ...
    "weight x4095";"raw gradient x64";"raw gradient x64";"R3"; ...
    "raw gradient x64";"weight codes per Q6 code x2^s";"shift bits";"bits";"Y full scale"];
values = table(Name,Value,Format,Unit,DeviceFormat);
end

function label = passLabel(passed)
if passed
    label = 'PASS';
else
    label = 'FAIL';
end
end
