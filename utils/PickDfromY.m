function D = PickDfromY(Y, Y_range, k)
	if nargin == 2 
		Y_range = [0, size(Y, 2)];
	end 
	C = numel(Y_range) - 1;
	D = [];
	for i = 1: C
		range = Y_range(i) + 1 : Y_range(i+1);
		Yi = Y(:, range);
		Ni = size(Yi,2);
		ids = randperm(Ni);
		D = [D, Yi(:, ids(1:k))];
	end
end