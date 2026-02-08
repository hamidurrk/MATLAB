function [minVal, maxVal, isMinFirst] = ex4p1(vec)
    minVal = min(vec);
    maxVal = max(vec);
    idxMin = find(vec == minVal, 1, 'first');
    idxMax = find(vec == maxVal, 1, 'first');
    isMinFirst = idxMin < idxMax;
end
