function W12 = divide49Round_fx(S49)
%DIVIDE49ROUND_FX 7×7 窗口和近似舍入除49，S49 为 [0,200655] 的 uint32。
% 49*2675=2^17+3，输出等于精确舍入均值或比它高 1 码。
x = uint32(S49) + uint32(24);                        % u18 [24,200679]
product = x .* uint32(2675);                         % u29 <=536816325
W12 = uint16(bitshift(product,-17));                 % u12 [0,4095]
end
