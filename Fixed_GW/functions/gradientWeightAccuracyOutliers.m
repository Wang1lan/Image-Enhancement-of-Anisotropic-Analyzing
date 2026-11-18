function rows = gradientWeightAccuracyOutliers(fx, ref, p)
%GRADIENTWEIGHTACCURACYOUTLIERS 超出整图验收限值的坐标、分支和误差来源。
names = ["Gx","Gy","G","W0","W7"];
limits = [0 0 0 1 2];
errors = {double(fx.Gx)-round(32760*ref.Gx), double(fx.Gy)-round(32760*ref.Gy), ...
    double(fx.G)-round(32760*ref.G), double(fx.W0_12)-round(4095*ref.W0), ...
    double(fx.W12)-round(4095*ref.W_smooth)};
rows = table;
for k = 1:numel(names)
    ix = find(abs(errors{k}) > limits(k));
    if isempty(ix) && width(rows)>0
        continue;
    end
    [PixelRow,PixelColumn] = ind2sub(size(fx.G),ix);
    Stage = repmat(names(k),numel(ix),1);
    ErrorLSB = errors{k}(ix);
    G = fx.G(ix);
    G_Q6 = fx.G_Q6(ix);
    Branch = repmat("middle",numel(ix),1);
    Branch(G_Q6<=p.Tlow_Q6) = "low";
    Branch(G_Q6>=p.Thigh_Q6) = "high";
    mid = Branch=="middle";
    Delta = nan(numel(ix),1);
    Delta(mid) = double(G_Q6(mid))-double(p.Tlow_Q6);
    C = repmat(double(p.K_Q17),numel(ix),1);
    s = repmat(17,numel(ix),1);
    W0FloatCodes = 4095*ref.W0(ix);
    W0Fixed = fx.W0_12(ix);
    S49 = fx.S49(ix);
    W7Fixed = fx.W12(ix);
    W7FloatCodes = 4095*ref.W_smooth(ix);
    W7CodeError = double(W7Fixed)-round(W7FloatCodes);
    w = double(p.Wmin12);
    endpoints = w+(4095-w)*min(max((double(G_Q6)-double(p.Tlow_Q6))/double(p.D_Q6),0),1);
    slopeCurve = endpoints;
    slopeCurve(mid) = min(w+Delta(mid)*double(p.K_Q17)/2^17,4095);
    EndpointError = endpoints-W0FloatCodes;
    SlopeError = slopeCurve-endpoints;
    PixelRoundError = double(W0Fixed)-slopeCurve;
    rows = [rows; table(Stage,PixelRow,PixelColumn,ErrorLSB,G,G_Q6,Branch,Delta,C,s, ...
        W0FloatCodes,W0Fixed,S49,W7Fixed,W7FloatCodes,W7CodeError, ...
        EndpointError,SlopeError,PixelRoundError)]; %#ok<AGROW>
end
end
