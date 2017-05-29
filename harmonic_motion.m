function W = harmonic_motion(X,PHASE,M)
  W = zeros(length(X),M);
  for t=[0:M-1]
    W(:,t+1) = X .* exp(-1i * 2 * pi * (t/(2*M) + PHASE));
  end
end