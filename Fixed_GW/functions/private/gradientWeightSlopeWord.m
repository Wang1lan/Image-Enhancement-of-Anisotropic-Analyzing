function K_Q17 = gradientWeightSlopeWord(Wmin12, D_Q6)
% 从已量化端点计算 R3，整组装载与离线标定共用 Q17 舍入规则。
span = uint64(4095) - uint64(Wmin12);                % u12 [2047,3685]
distance = uint64(D_Q6);                             % u18 [6289,209664]
halfD = bitshift(distance, -1);                      % u17 [3144,104832]
numerator = bitshift(span,17) + halfD;                % u29 <=483105152
K_Q17 = uint32(idivide(numerator,distance,'floor'));  % u17 [1280,76801]
end
