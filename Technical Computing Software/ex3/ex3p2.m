function T = ex3p2(C)
	validateattributes(C, {'numeric'}, {'vector'});
	C = C(:); % ensure column vector

	F = (9/5) .* C + 32;

	T = table(C, F, 'VariableNames', {'Celsius', 'Fahrenheit'});
end

