% A useful function implementing the jackknife error formula.
% 
% Each column of the input represents a different measured quantity. 
% Each row of the input represents a different jackknife sample
% The first row of the input should be the entire sample, from which
% the central value is computed

function [central_values, jackknife_errors] = JackknifeError(jackknife_values)
  central_values = jackknife_values(1, :);
  if size(jackknife_values, 1) == 1
    jackknife_errors = zeros(1, size(jackknife_values, 2));
    return
  end

  subsample_values = jackknife_values(2:end, :);
  N = size(jackknife_values, 1)-1;
  jackknife_errors = sqrt(N/(N-1) * sum((subsample_values - repmat(central_values, N, 1)).^2, 1));
end
