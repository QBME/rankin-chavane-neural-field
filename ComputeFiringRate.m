function [S,DS] = ComputeFiringRate(v,th)

  S = 1./ (1 + exp( - v + th) ) - 1./ (1 + exp(th) );

  if nargout > 1
    DS = exp( - v + th)./ (1 + exp( - v + th) ).^2;
  end

end
