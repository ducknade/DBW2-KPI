function DoFitJackknife

  %TODO: generalize to non-16x32

  % Compile C++ functions
  mex ComputePeriodicModelResiduals.C
  mex ComputeOpenModelResiduals.C

  %%% Periodic ensembles %%%
%  guess_params = [18, 0.115];
%  Njack = 97;
%  [central_params, param_errors, chisq, chisq_error] = RunEnsemble('periodic-8x16', 16, 'periodic', Njack, guess_params);
%  WritePeriodicResults('periodic-8x16', 16, central_params, param_errors, chisq, chisq_error);
%
%  guess_params = [18, 0.115];
%  Njack = 265;
%  [central_params, param_errors, chisq, chisq_error] = RunEnsemble('shorttraj-periodic-8x16', 16, 'periodic', Njack, guess_params);
%  WritePeriodicResults('shorttraj-periodic-8x16', 16, central_params, param_errors, chisq, chisq_error);
%
%  guess_params = [18, 0.115];
%  Njack = 50;
%  [central_params, param_errors, chisq, chisq_error] = RunEnsemble('longtraj-periodic-8x16', 16, 'periodic', Njack, guess_params);
%  WritePeriodicResults('longtraj-periodic-8x16', 16, central_params, param_errors, chisq, chisq_error);
% 
%  guess_params = [18, 0.115];
%  Njack = 180;
%  [central_params, param_errors, chisq, chisq_error] = RunEnsemble('longtraj2-periodic-8x16', 16, 'periodic', Njack, guess_params);
%  WritePeriodicResults('longtraj2-periodic-8x16', 16, central_params, param_errors, chisq, chisq_error);
%
%  guess_params = [55, 0.115];
%  Njack = 122;
%  [central_params, param_errors, chisq, chisq_error] = RunEnsemble('periodic-10x20', 20, 'periodic', Njack, guess_params);
%  WritePeriodicResults('periodic-10x20', 20, central_params, param_errors, chisq, chisq_error);
%
%  guess_params = [200, 0.115];
%  Njack = 54;
%  [central_params, param_errors, chisq, chisq_error] = RunEnsemble('periodic-12x24', 24, 'periodic', Njack, guess_params);
%  WritePeriodicResults('periodic-12x24', 24, central_params, param_errors, chisq, chisq_error);
%
%  guess_params = [550, 0.115];
%  Njack = 37;
%  [central_params, param_errors, chisq, chisq_error] = RunEnsemble('periodic-14x28', 28, 'periodic', Njack, guess_params);
%  WritePeriodicResults('periodic-14x28', 28, central_params, param_errors, chisq, chisq_error);
%
%  guess_params = [2200, 0.115];
%  Njack = 21;
%  [central_params, param_errors, chisq, chisq_error] = RunEnsemble('periodic-16x32', 32, 'periodic', Njack, guess_params);
%  WritePeriodicResults('periodic-16x32', 32, central_params, param_errors, chisq, chisq_error);
%
%  guess_params = [2200, 0.115];
%  Njack = 17;
%  [central_params, param_errors, chisq, chisq_error] = RunEnsemble('iwasaki-closed-2', 64, 'periodic', Njack, guess_params);
%  WritePeriodicResults('iwasaki-closed-2', 64, central_params, param_errors, chisq, chisq_error);


  %%% Open ensembles %%%
  guess_params = [18, repmat(0.115, 1, 8)];
  Njack = 88;
  [central_params, param_errors, chisq, chisq_error] = RunEnsemble('open-8x16', 16, 'open', Njack, guess_params);
  WriteOpenResults('open-8x16', 16, central_params, param_errors, chisq, chisq_error);

  guess_params = [55, repmat(0.115, 1, 10)];
  Njack = 121;
  [central_params, param_errors, chisq, chisq_error] = RunEnsemble('open-10x20', 20, 'open', Njack, guess_params);
  WriteOpenResults('open-10x20', 20, central_params, param_errors, chisq, chisq_error);

  guess_params = [200, repmat(0.115, 1, 12)];
  Njack = 51;
  [central_params, param_errors, chisq, chisq_error] = RunEnsemble('open-12x24', 24, 'open', Njack, guess_params);
  WriteOpenResults('open-12x24', 24, central_params, param_errors, chisq, chisq_error);

  guess_params = [550, repmat(0.115, 1, 14)];
  Njack = 49;
  [central_params, param_errors, chisq, chisq_error] = RunEnsemble('open-14x28', 28, 'open', Njack, guess_params);
  WriteOpenResults('open-14x28', 28, central_params, param_errors, chisq, chisq_error);

  guess_params = [2300, repmat(0.115, 1, 16)];
  Njack = 22;
  [central_params, param_errors, chisq, chisq_error] = RunEnsemble('open-16x32', 32, 'open', Njack, guess_params);
  WriteOpenResults('open-16x32', 32, central_params, param_errors, chisq, chisq_error);
end


function [central_params, param_errors, chisq, chisq_error] = RunEnsemble(ens_name, T, bc, Njack, guess_params)
  fprintf('Fitting ensemble %s; %s boundary conditions; %d jackknife blocks\n', ens_name, bc, Njack);

  if strcmp(bc, 'open')
    model_residuals = @ComputeOpenModelResiduals;
  elseif strcmp(bc, 'periodic')
    model_residuals = @ComputePeriodicModelResiduals;
  else 
    error('Bad bc input: %s', bc);
  end

  jackknife_params = [];
  jackknife_chisq = [];
  for jack = 0:Njack-1
    fprintf('Fitting jackknife #%d\n', jack);
    [jackknife_params(jack+1,:), jackknife_chisq(jack+1,1), res] = ...
        lsqnonlin(@(x) model_residuals(ens_name, T, jack, x), guess_params, [], [], optimset('Display', 'none'));
  end

  % Do shape:
  tmp = model_residuals(ens_name, T, 0, jackknife_params(1,:), 'shape');

  [central_params, param_errors] = JackknifeError(jackknife_params);
  [chisq, chisq_error] = JackknifeError(jackknife_chisq);
end

function WritePeriodicResults(ens_name, T, central_params, param_errors, chisq, chisq_error)
  dlmwrite(sprintf('results/%s-tauexp.dat', ens_name), [central_params(1), param_errors(1), T], '\t');
  dlmwrite(sprintf('results/%s-D.dat', ens_name), [central_params(2), param_errors(2), T], '\t');

  dlmwrite(sprintf('results/%s-chisq.dat', ens_name), [chisq, chisq_error], '\t');
end

function WriteOpenResults(ens_name, T, central_params, param_errors, chisq, chisq_error)
  dlmwrite(sprintf('results/%s-tauexp.dat', ens_name), [central_params(1), param_errors(1), T], '\t');
%  data = [(0:T/2-1)', central_params(2:end)', param_errors(2:end)', repmat(T, length(central_params)-1, 1)];
  data = [0.5+(0:T/2-1)', central_params(2:end)', param_errors(2:end)', repmat(T, length(central_params)-1, 1)];
  dlmwrite(sprintf('results/%s-D.dat', ens_name), data, '\t');

  dlmwrite(sprintf('results/%s-chisq.dat', ens_name), [chisq, chisq_error], '\t');
end

