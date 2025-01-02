function [Nnodes,Nel,NDOFs] = GetDiscretization(xn,Tn)
% Some variables
Nnodes = length(xn);
Nel = size(Tn,1);
NDOFs = Nnodes*6;

end

