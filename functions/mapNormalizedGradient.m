function y = mapNormalizedGradient(t, mapping)
%MAPNORMALIZEDGRADIENT 将阈值归一化梯度映射为 [0,1] 系数。
%   y = mapNormalizedGradient(t,mapping)
%   t 已裁剪到 [0,1]；mapping 已由 computeGradientWeight 校验。
%   输出与 t 同尺寸；gamma 仅在 power 模式参与计算。
switch mapping.type
    case 'linear'
        y = t;
    case 'power'
        y = t .^ mapping.gamma;
    case 'smoothstep'
        y = t .* t .* (3 - 2 .* t);
    otherwise
        error('mapNormalizedGradient:MappingType', ...
            'Unsupported mapping type: %s', mapping.type);
end
end
