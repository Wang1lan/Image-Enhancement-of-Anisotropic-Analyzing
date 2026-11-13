function [W0_12, G_Q4] = mapGradientWeight_fx(G, params)
%MAPGRADIENTWEIGHT_FX 原始整数 Sobel L1 幅值映射为 u12 初始权重。
%   [W0_12,G_Q4] = mapGradientWeight_fx(G,params)
%   G 为 uint16 容器中的 u15，范围 [0,24570]；params 由预处理函数生成。
%   W0_12 与 G 同尺寸，uint16，范围 [params.Wmin12,4095]。

G_Q4 = bitshift(uint32(G), 4);                       % u19 [0,393120]
W0_12 = repmat(params.Wmin12, size(G));               % u12 [Wmin12,4095]
W0_12(G_Q4 >= params.Thigh_Q4) = uint16(4095);
middle = G_Q4 > params.Tlow_Q4 & G_Q4 < params.Thigh_Q4;

% 仅中间区间参与减法，故 delta 非负且严格小于 D_Q4。
delta_Q4 = uint64(G_Q4(middle)) - uint64(params.Tlow_Q4); % u19 [0,D_Q4-1]
product_u33 = delta_Q4 .* uint64(params.K_Q20);       % u33 规格预留，<=4294180800
roundedNum_u33 = product_u33 + uint64(524288);        % u33 规格预留，<=4294705088
increment_u12 = bitshift(roundedNum_u33, -20);         % u12 [0,4095-Wmin12]
weightSum_u13 = uint64(params.Wmin12) + increment_u12; % u13 规格预留，保守 [0,8190]
weightClip_u12 = min(max(weightSum_u13, uint64(params.Wmin12)), ...
    uint64(4095));                                  % u12 [Wmin12,4095]
W0_12(middle) = uint16(weightClip_u12);                % u12，显式限幅后缩窄
end
