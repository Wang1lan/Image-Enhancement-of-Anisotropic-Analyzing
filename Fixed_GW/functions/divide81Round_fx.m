function W12 = divide81Round_fx(S81)
%DIVIDE81ROUND_FX 对 u19 九乘九累加和执行精确的最近整数除81。
%   W12 = divide81Round_fx(S81)
%   S81 为 uint32，范围 [0,331695]；输出 uint16，范围 [0,4095]。
%   全合法域上 ((S81+40)*828505)>>26 等于 floor((S81+40)/81)。

sumRound_u19 = uint64(S81) + uint64(40);             % u19 [40,331735]
meanNum_u38 = sumRound_u19 .* uint64(828505);        % u38 [33140200,274844106175]
mean_u12 = bitshift(meanNum_u38, -26);                % u12 [0,4095]
meanClip_u12 = min(max(mean_u12, uint64(0)), uint64(4095)); % u12
W12 = uint16(meanClip_u12);                          % u12，显式限幅后缩窄
end
