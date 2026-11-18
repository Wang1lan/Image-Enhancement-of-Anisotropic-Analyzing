function [W12, debugFx] = computeGradientWeight_fx(Y12, params)
%COMPUTEGRADIENTWEIGHT_FX 12-bit Y -> 12-bit 梯度自适应权重，纯整数。
%   [W12,debugFx] = computeGradientWeight_fx(Y12,params)
%   Y12 为至少 3x3 的二维 uint16，范围 [0,4095]。
%   params 由 prepareGradientWeightFxParams 或 loadGradientWeightFxRegs 生成；输出同尺寸 uint16。
%   Sobel 及 7x7 均值均采用 replicate，Sobel 符号遵从仓库卷积定义。
%   第二输出按需返回整数阶段矩阵；uN/sN 均为总位宽，sN 已包含符号位。

assert(isa(Y12, 'uint16') && isreal(Y12) && ismatrix(Y12) && ...
    all(size(Y12) >= 3) && all(Y12(:) <= uint16(4095)), ...
    'computeGradientWeight_fx:Input', ...
    'Y12 必须是至少 3x3、数值不超过4095的二维 uint16 图像。');
[height, width] = size(Y12);
pixels_u12 = int32(padarray(Y12, [1 1], 'replicate', 'both')); % u12 输入

% 翻转 Kx 后卷积，左减右。差值可为负，所有项先提升到 int32。
xTop_s13 = pixels_u12(1:height, 1:width) - ...
    pixels_u12(1:height, 3:width+2);                  % s13 [-4095,4095]
xMiddle_s13 = pixels_u12(2:height+1, 1:width) - ...
    pixels_u12(2:height+1, 3:width+2);                % s13 [-4095,4095]
xBottom_s13 = pixels_u12(3:height+2, 1:width) - ...
    pixels_u12(3:height+2, 3:width+2);                % s13 [-4095,4095]
xMiddle_s14 = xMiddle_s13 .* int32(2);               % s14 [-8190,8190]
xPartial_s15 = xTop_s13 + xMiddle_s14;               % s15 [-12285,12285]
Gx = xPartial_s15 + xBottom_s13;                     % s15 [-16380,16380]

yLeft_s13 = pixels_u12(1:height, 1:width) - ...
    pixels_u12(3:height+2, 1:width);                  % s13 [-4095,4095]
yMiddle_s13 = pixels_u12(1:height, 2:width+1) - ...
    pixels_u12(3:height+2, 2:width+1);                % s13 [-4095,4095]
yRight_s13 = pixels_u12(1:height, 3:width+2) - ...
    pixels_u12(3:height+2, 3:width+2);                % s13 [-4095,4095]
yMiddle_s14 = yMiddle_s13 .* int32(2);               % s14 [-8190,8190]
yPartial_s15 = yLeft_s13 + yMiddle_s14;               % s15 [-12285,12285]
Gy = yPartial_s15 + yRight_s13;                      % s15 [-16380,16380]

xAbs_u14 = uint32(abs(Gx));                         % u14 [0,16380]，int32 取绝对值安全
yAbs_u14 = uint32(abs(Gy));                         % u14 [0,16380]
gradient_u15 = xAbs_u14 + yAbs_u14;                  % u15 [0,24570]，共享邻域紧上界
gradientClip_u15 = min(gradient_u15, uint32(24570));   % u15 [0,24570]
G = uint16(gradientClip_u15);                         % u15，范围证明及限幅后缩窄
[W0_12, G_Q6] = mapGradientWeight_fx(G, params);

% 两个方向的累加均保留全精度，到49项齐全后只舍入一次。
horizontalPad = padarray(uint32(W0_12), [0 3], 'replicate', 'both');
S7 = zeros(height, width, 'uint32');                  % u15 [0,28665]
for offset = 1:7
    S7 = S7 + horizontalPad(:, offset:offset+width-1);
end
verticalPad = padarray(S7, [3 0], 'replicate', 'both');
S49 = zeros(height, width, 'uint32');                 % u18 [0,200655]
for offset = 1:7
    S49 = S49 + verticalPad(offset:offset+height-1, :);
end
W12 = divide49Round_fx(S49);                          % u12 [Wmin12,4095]
if nargout > 1
    debugFx = struct('Gx', Gx, 'Gy', Gy, 'G', G, 'G_Q6', G_Q6, ...
        'W0_12', W0_12, 'S7', S7, 'S49', S49, 'W12', W12);
end
end
