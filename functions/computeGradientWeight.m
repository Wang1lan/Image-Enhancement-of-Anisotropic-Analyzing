function [W, debug] = computeGradientWeight(I, opts)
%COMPUTEGRADIENTWEIGHT 保边降噪引导的 Sobel L1 梯度增强权重。
%   [W,debug] = computeGradientWeight(I,opts)
%   I：至少 3x3 的二维 uint8/uint16/single/double，范围 [0,inputMax]。
%   GLow/GHigh 必填，单位为固定满量程归一化后的 Sobel L1 幅值。
%   有效区 W 属于 [Wmin,1]；validMask 只屏蔽最终 W 和统计，不改变滤波输入。
%   denoise.sigmaRange/guidedEps 使用输入原生灰度级/灰度级平方。
%   debugLevel='off'/'stats'/'maps'；off 或不取第二输出时跳过附加统计。
%   本函数只分析亮度，不修改增强底图或带通细节，不绘图或保存文件。

p = parseOptions(opts);
validateattributes(I, {'uint8','uint16','single','double'}, ...
    {'real','2d','finite','nonnegative','<=',p.inputMax}, mfilename, 'I');
assert(all(size(I) >= 3), 'computeGradientWeight:InputSize', ...
    '亮度图 I 必须至少为 3x3。');
Yinput = double(I) / p.inputMax;
if isempty(p.validMask)
    mask = true(size(I));
else
    assert(islogical(p.validMask) && isequal(size(p.validMask), size(I)), ...
        'computeGradientWeight:Mask', 'validMask 必须为空或与 I 同尺寸的 logical。');
    mask = p.validMask;
end

Yguide = Yinput;
guideBoundary = 'bypass';
if p.denoise.enable && p.denoise.windowSize > 1
    d = p.denoise;
    switch d.method
        case 'bilateral'
            Yguide = bilateralFilterShiftSum(Yinput, d.windowSize, ...
                d.sigmaSpatial, d.sigmaRange / p.inputMax);
            guideBoundary = 'replicate';
        case 'guided'
            assert(d.windowSize <= min(size(I)), ...
                'computeGradientWeight:GuidedWindow', '引导滤波窗口不能超过图像尺寸。');
            assert(exist('imguidedfilter', 'file') == 2, ...
                'computeGradientWeight:GuidedDependency', '引导降噪需要 imguidedfilter。');
            Yguide = imguidedfilter(Yinput, Yinput, ...
                'NeighborhoodSize', [d.windowSize d.windowSize], ...
                'DegreeOfSmoothing', d.guidedEps / p.inputMax^2);
            guideBoundary = 'MATLAB_builtin';
    end
end
assert(all(isfinite(Yguide(:))), 'computeGradientWeight:GuideFinite', ...
    '引导降噪结果含非有限值。');
Yguide = min(max(Yguide, 0), 1);

% filterShiftSum 执行卷积；梯度符号遵从卷积定义，L1 幅值不受符号影响。
Kx = [-1 0 1; -2 0 2; -1 0 1] / 8;
Gx = filterShiftSum(Yguide, Kx, 'replicate');
Gy = filterShiftSum(Yguide, Kx', 'replicate');
G = abs(Gx) + abs(Gy);
t = min(max((G-p.GLow) / (p.GHigh-p.GLow), 0), 1);
W = p.Wmin + (1-p.Wmin) .* t;
W(~mask) = 0;

debug = struct();
if nargout < 2 || strcmp(p.debugLevel, 'off')
    return;
end
debug.paramsUsed = rmfield(p, 'validMask');
debug.paramsUsed.sobelScale = 1/8;
debug.paramsUsed.gradientDefinition = 'abs(Gx)+abs(Gy)';
debug.paramsUsed.gradientConvention = 'convolution';
debug.paramsUsed.gradientBoundary = 'replicate';
debug.paramsUsed.guideBoundary = guideBoundary;
debug.messages = {};
debug.stats.validCount = nnz(mask);
debug.stats.totalCount = numel(mask);
debug.stats.G = distributionStats(G(mask));
debug.stats.W = distributionStats(W(mask));
debug.stats.ratios = struct('GLow', mean(G(mask) <= p.GLow), ...
    'GHigh', mean(G(mask) >= p.GHigh), 'WMin', mean(W(mask) == p.Wmin), ...
    'WHigh', mean(W(mask) >= .8));
if ~any(mask(:))
    debug.messages{end+1} = 'No valid pixels: count=0; statistics and ratios are NaN.';
end
if strcmp(p.debugLevel, 'maps')
    debug.Yinput = Yinput;
    debug.Yguide = Yguide;
    debug.Gx = Gx;
    debug.Gy = Gy;
    debug.G = G;
    debug.W = W;
    debug.validMask = mask;
end
end

function p = parseOptions(opts)
denoiseDefaults = struct('enable', true, 'method', 'bilateral', ...
    'windowSize', 5, 'sigmaSpatial', 1.5, 'sigmaRange', 10, 'guidedEps', 4);
p = struct('inputMax', 255, 'GLow', [], 'GHigh', [], 'Wmin', .3, ...
    'validMask', [], 'debugLevel', 'off', 'denoise', denoiseDefaults);
p = mergeOptions(p, opts);
validateattributes(p.inputMax, {'numeric'}, {'real','scalar','finite','positive'});
validateattributes(p.Wmin, {'numeric'}, {'real','scalar','finite','>=',0,'<=',1});
p.inputMax = double(p.inputMax);
p.Wmin = double(p.Wmin);
thresholds = {p.GLow, p.GHigh};
for k = 1:2
    value = thresholds{k};
    assert(isnumeric(value) && isreal(value) && isscalar(value) && isfinite(value), ...
        'computeGradientWeight:Thresholds', ...
        'GLow/GHigh 必须为有限实数标量；请先运行 demo_gradientWeight，再填写标定门限。');
end
assert(0 <= p.GLow && p.GLow < p.GHigh && p.GHigh <= 1, ...
    'computeGradientWeight:Thresholds', '门限要求 0 <= GLow < GHigh <= 1。');
p.GLow = double(p.GLow);
p.GHigh = double(p.GHigh);
p.debugLevel = validatestring(p.debugLevel, {'off','stats','maps'});
p.denoise = mergeOptions(denoiseDefaults, p.denoise);
d = p.denoise;
validateattributes(d.enable, {'logical'}, {'scalar'});
if d.enable
    d.method = validatestring(d.method, {'bilateral','guided'});
    validateattributes(d.windowSize, {'numeric'}, ...
        {'real','scalar','finite','integer','positive','odd'});
    d.windowSize = double(d.windowSize);
    if strcmp(d.method, 'bilateral')
        validateattributes(d.sigmaSpatial, {'numeric'}, {'real','scalar','finite','positive'});
        validateattributes(d.sigmaRange, {'numeric'}, {'real','scalar','finite','positive'});
        d.sigmaSpatial = double(d.sigmaSpatial);
        d.sigmaRange = double(d.sigmaRange);
    else
        validateattributes(d.guidedEps, {'numeric'}, {'real','scalar','finite','positive'});
        d.guidedEps = double(d.guidedEps);
    end
end
p.denoise = d;
end

function p = mergeOptions(p, opts)
assert(isstruct(opts) && isscalar(opts), ...
    'computeGradientWeight:Options', 'opts 和 denoise 必须为标量结构体。');
names = fieldnames(opts);
unknown = setdiff(names, fieldnames(p));
assert(isempty(unknown), 'computeGradientWeight:UnknownOption', '存在未知选项字段。');
for k = 1:numel(names)
    p.(names{k}) = opts.(names{k});
end
end

function st = distributionStats(v)
n = numel(v);
st = struct('count', n, 'min', NaN, 'max', NaN, 'mean', NaN, 'std', NaN, ...
    'P05', NaN, 'P50', NaN, 'P75', NaN, 'P90', NaN, 'P95', NaN, 'P99', NaN);
if n == 0
    return;
end
v = sort(v(:)).';
st.min = v(1);
st.max = v(end);
st.mean = mean(v);
st.std = std(v, 1);
pos = 1+(n-1)*[5 50 75 90 95 99]/100;
q = v(floor(pos)) + (pos-floor(pos)).*(v(ceil(pos))-v(floor(pos)));
names = {'P05','P50','P75','P90','P95','P99'};
for k = 1:numel(names)
    st.(names{k}) = q(k);
end
end
