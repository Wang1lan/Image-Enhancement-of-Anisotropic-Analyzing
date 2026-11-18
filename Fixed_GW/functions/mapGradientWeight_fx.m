function [W0_12, G_Q6] = mapGradientWeight_fx(G, params)
%MAPGRADIENTWEIGHT_FX 将原始 Sobel L1 梯度映射为 u12 初始权重。
%   G 为 [0,24570] 的 uint16；params 为单帧完整配置。G_Q6 保留全部 u21。
G_Q6 = bitshift(uint32(G),6);                         % u21 [0,1572480]
W0_12 = repmat(params.Wmin12,size(G));                 % u12
W0_12(G_Q6 >= params.Thigh_Q6) = uint16(4095);
middle = G_Q6 > params.Tlow_Q6 & G_Q6 < params.Thigh_Q6;
delta = G_Q6(middle) - params.Tlow_Q6;                 % u18 [1,D_Q6-1]
% 整组寄存器已校验，以下乘积与偏置之和 <=483170688，uint32 足以承载 u29。
product = delta .* params.K_Q17;                      % u29 <=483105152
roundedNum = product + uint32(65536);                 % u29 <=483170688，偏置 u17
increment = bitshift(roundedNum,-17);                 % u12 <=3686
span = uint32(4095) - uint32(params.Wmin12);          % u12
W0_12(middle) = uint16(uint32(params.Wmin12) + min(increment,span)); % u12
end
