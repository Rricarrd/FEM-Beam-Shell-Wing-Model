function [Nnodes,Nel,NDOFs] = GetDiscretization(xn)
% Some variables
Nnodes = length(xn);
Nel = Nnodes - 1;
NDOFs = Nnodes*6;
end

