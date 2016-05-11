function A = myForm(M, N, k)
	I = eye(k);
	A = kron(I, M) + repmat(N, k, k);
end 