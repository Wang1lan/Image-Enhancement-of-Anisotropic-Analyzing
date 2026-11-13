function params = prepareGradientWeightFxParams(GLow, GHigh, Wmin)
%PREPAREGRADIENTWEIGHTFXPARAMS 三个归一化参数转换为整数权重寄存器。
%   params = prepareGradientWeightFxParams(GLow,GHigh,Wmin)
%   0 <= GLow < GHigh <= 1，0 <= Wmin <= 1，且 Q4 量化后门限不同。
%   本函数为上位机配置阶段，允许浮点运算及整数除法。
%   uN/sN 均指总位宽，sN 已包含符号位。

values = {GLow, GHigh, Wmin};
for index = 1:numel(values)
    value = values{index};
    assert(isnumeric(value) && isreal(value) && isscalar(value) && ...
        isfinite(value) && value >= 0 && value <= 1, ...
        'prepareGradientWeightFxParams:Parameters', ...
        'GLow、GHigh、Wmin 必须是 [0,1] 内的有限实数标量。');
end
assert(GLow < GHigh, 'prepareGradientWeightFxParams:Parameters', ...
    '门限必须满足 GLow < GHigh。');

% 归一化 Sobel L1 -> 原始 Sobel Q4，524160 = 8*4095*16。
thresholds = round([double(GLow), double(GHigh)] * 524160);
assert(thresholds(1) < thresholds(2), ...
    'prepareGradientWeightFxParams:Parameters', 'Q4 量化后的高低门限必须不同。');
Tlow_Q4 = uint32(min(max(thresholds(1), 0), 524160)); % u19 [0,524160]
Thigh_Q4 = uint32(min(max(thresholds(2), 0), 524160)); % u19 [0,524160]
D_Q4 = uint32(uint64(Thigh_Q4) - uint64(Tlow_Q4));    % u19 [1,524160]
weightRounded = round(double(Wmin) * 4095);           % 配置阶段最近舍入，半值向上
Wmin12 = uint16(min(max(weightRounded, 0), 4095));     % u12 [0,4095]

weightSpan_u12 = uint64(4095) - uint64(Wmin12);       % u12 [0,4095]
slopeNum_u32 = bitshift(weightSpan_u12, 20);          % u32 [0,4293918720]
halfD_u18 = idivide(uint64(D_Q4), uint64(2), 'floor'); % u18 [0,262080]
slopeRound_u32 = slopeNum_u32 + halfD_u18;            % u32 [0,4294180800]
slope_u32 = idivide(slopeRound_u32, uint64(D_Q4), 'floor'); % u32 [0,4293918720]
assert(slope_u32 <= uint64(4294967295), ...
    'prepareGradientWeightFxParams:Parameters', 'K_Q20 超出 u32 寄存器范围。');
K_Q20 = uint32(min(slope_u32, uint64(4294967295)));    % u32

params = struct('GLow', double(GLow), 'GHigh', double(GHigh), ...
    'Wmin', double(Wmin), 'Tlow_Q4', Tlow_Q4, 'Thigh_Q4', Thigh_Q4, ...
    'D_Q4', D_Q4, 'Wmin12', Wmin12, 'K_Q20', K_Q20, ...
    'fracBitsGradient', 4, 'fracBitsSlope', 20, 'inputMax', 4095);
end
