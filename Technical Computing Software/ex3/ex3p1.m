function [sumVec, diffVec, prodVec, divVec] = ex3p1(x)
	validateattributes(x, {'numeric'}, {'vector', 'numel', 5});
	x = reshape(x,1,5);
	y = 1:5; % integer array on closed interval [1,5]

	sumVec  = x + y;
	diffVec = x - y;
	prodVec = x .* y;
	divVec  = x ./ y;
end

