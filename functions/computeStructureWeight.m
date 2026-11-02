function [W, debug] = computeStructureWeight(I, opts)
%COMPUTESTRUCTUREWEIGHT 由方向一致性与梯度能量生成 [0,1] 增强权重。
%   W = computeStructureWeight(I, opts)
%   [W, debug] = computeStructureWeight(I, opts)
%   I：至少 3x3 的二维 uint8/uint16/single/double，实数、有限、非负。
%   固定标度 Y=double(I)/inputMax；不得用逐图最大值替代 inputMax。
%   opts.SLow、SHigh 必填，单位为归一化输入的中心差分梯度平方。
%   其他默认值：inputMax=1，sigma=1，rho=2，CLow=.2，CHigh=.8，
%   epsilon=1e-12，validMask=[]（全有效），debugLevel='off'。
%   validMask 只控制最终 W 和统计取样，不改变滤波输入。
%   debugLevel：off 返回空结构体；stats 返回 paramsUsed/stats/messages；
%   maps 另返回分析图。未请求第二输出时跳过 debug 构建与附加统计。
%   stats.C/S/WC/WS/W 含 count/min/max/mean/std/P05/P50/P75/P90/P95/P99；
%   ratios.CLow/CHigh/SLow/SHigh 分别为 <=低门限、>=高门限的比例，
%   ratios.WZero/WHigh 分别为 W==0、W>=.8 的比例，仅统计有效像素。
%   本函数只计算权重；增强由调用方执行 base + gain.*W.*detail。

p = parseOptions(opts);
if ~(isa(I, 'uint8') || isa(I, 'uint16') || isa(I, 'single') || isa(I, 'double')) ...
        || ~isreal(I) || ~ismatrix(I) || any(size(I) < 3)
    error('computeStructureWeight:Input', ...
        'I must be a real 2-D uint8/uint16/single/double image, at least 3x3.');
end
Y = double(I);
if any(~isfinite(Y(:)) | Y(:) < 0 | Y(:) > p.inputMax)
    error('computeStructureWeight:InputRange', 'I must be finite and within [0,inputMax].');
end
if isempty(p.validMask)
    mask = true(size(Y));
elseif ~islogical(p.validMask) || ~isequal(size(p.validMask), size(Y))
    error('computeStructureWeight:Mask', 'validMask must be empty or a logical mask matching I.');
else
    mask = p.validMask;
end
Y = Y / p.inputMax;

Ysmooth = gaussianReplicate(Y, p.sigma);
Ix = (Ysmooth(:, [2:end end]) - Ysmooth(:, [1 1:end-1])) / 2;
Iy = (Ysmooth([2:end end], :) - Ysmooth([1 1:end-1], :)) / 2;
Jxx = gaussianReplicate(Ix.^2, p.rho);
Jyy = gaussianReplicate(Iy.^2, p.rho);
Jxy = gaussianReplicate(Ix.*Iy, p.rho);
S = max(Jxx + Jyy, 0);
R = sqrt(max((Jxx-Jyy).^2 + 4*Jxy.^2, 0));
C = min(max(R ./ (S+p.epsilon), 0), 1);
C(S == 0) = 0;
WC = min(max((C-p.CLow) / (p.CHigh-p.CLow), 0), 1);
WS = min(max((S-p.SLow) / (p.SHigh-p.SLow), 0), 1);
W = WC .* WS;
W(~mask) = 0;

debug = struct();
if nargout < 2 || strcmp(p.debugLevel, 'off')
    return;
end
debug.paramsUsed = rmfield(p, 'validMask');
debug.paramsUsed.gradientOperator = 'central_difference_half';
debug.paramsUsed.boundary = 'replicate';
debug.paramsUsed.gaussianRadius = 'ceil(3*scale)';
debug.messages = {};
if p.rho == 0
    debug.messages{end+1} = 'rho=0: single-point gradient outer products are rank one; C may approach 1.';
end
debug.stats.validCount = nnz(mask);
debug.stats.totalCount = numel(Y);
if ~any(mask(:))
    debug.messages{end+1} = 'No valid pixels: distribution statistics and ratios are NaN.';
end
debug.stats.C = distributionStats(C(mask));
debug.stats.S = distributionStats(S(mask));
debug.stats.WC = distributionStats(WC(mask));
debug.stats.WS = distributionStats(WS(mask));
debug.stats.W = distributionStats(W(mask));
debug.stats.ratios = struct('CLow', mean(C(mask) <= p.CLow), ...
    'CHigh', mean(C(mask) >= p.CHigh), 'SLow', mean(S(mask) <= p.SLow), ...
    'SHigh', mean(S(mask) >= p.SHigh), 'WZero', mean(W(mask) == 0), ...
    'WHigh', mean(W(mask) >= .8));
if strcmp(p.debugLevel, 'maps')
    debug.Y = Y;
    debug.Ysmooth = Ysmooth;
    debug.Ix = Ix;
    debug.Iy = Iy;
    debug.Jxx = Jxx;
    debug.Jyy = Jyy;
    debug.Jxy = Jxy;
    debug.C = C;
    debug.S = S;
    debug.WC = WC;
    debug.WS = WS;
    debug.validMask = mask;
    debug.W = W;
end
end

function p = parseOptions(opts)
p = struct('inputMax', 1, 'sigma', 1, 'rho', 2, 'CLow', .2, 'CHigh', .8, ...
    'SLow', [], 'SHigh', [], 'epsilon', 1e-12, 'validMask', [], 'debugLevel', 'off');
if ~isstruct(opts) || ~isscalar(opts)
    error('computeStructureWeight:Options', 'opts must be a scalar struct.');
end
names = fieldnames(opts);
unknown = setdiff(names, fieldnames(p));
if ~isempty(unknown)
    error('computeStructureWeight:UnknownOption', 'Unknown option: %s.', unknown{1});
end
for k = 1:numel(names)
    p.(names{k}) = opts.(names{k});
end
scalarNames = {'inputMax', 'sigma', 'rho', 'CLow', 'CHigh', 'SLow', 'SHigh', 'epsilon'};
for k = 1:numel(scalarNames)
    name = scalarNames{k};
    value = p.(name);
    if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ~isfinite(value)
        error('computeStructureWeight:Parameter', '%s must be a finite real scalar; SLow/SHigh are required.', name);
    end
    p.(name) = double(value);
end
if p.inputMax <= 0 || p.sigma < 0 || p.rho < 0 || p.epsilon <= 0
    error('computeStructureWeight:ParameterRange', 'inputMax/epsilon must be positive; sigma/rho must be nonnegative.');
end
if ~(0 <= p.CLow && p.CLow < p.CHigh && p.CHigh <= 1) ...
        || ~(0 <= p.SLow && p.SLow < p.SHigh)
    error('computeStructureWeight:Thresholds', 'Require 0<=CLow<CHigh<=1 and 0<=SLow<SHigh.');
end
level = p.debugLevel;
if ~((ischar(level) && isrow(level)) || (isstring(level) && isscalar(level))) ...
        || ~any(strcmp(level, {'off', 'stats', 'maps'}))
    error('computeStructureWeight:DebugLevel', 'debugLevel must be off, stats or maps.');
end
p.debugLevel = char(level);
end

function B = gaussianReplicate(A, scale)
if scale == 0
    B = A;
    return;
end
r = ceil(3*scale);
g = exp(-.5*((-r:r)/scale).^2);
g = g/sum(g);
cols = min(max((1-r):(size(A,2)+r), 1), size(A,2));
B = conv2(A(:, cols), g, 'valid');
rows = min(max((1-r):(size(A,1)+r), 1), size(A,1));
B = conv2(B(rows, :), g', 'valid');
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
names = {'P05', 'P50', 'P75', 'P90', 'P95', 'P99'};
for k = 1:numel(names)
    st.(names{k}) = q(k);
end
end
