function [loss, d_loss] = GSM_mvn_gradfit_loss(state, data, params)

model = [];
model.pis = exp(state.log_pis*20);
model.pis = model.pis ./ sum(model.pis);
model.Sigma = vec2PSD(state.Sigma_vec);
model.vars = data.vars;
if isfield(state, 'mu')
  model.mu = state.mu;
end

loss = -sum(GSM_mvn_pdf(model, data.X));

% try
%   P = 0;
%   for j = 1:length(model.vars)
%     P = P + exp(log(model.pis(j)) + lmvnpdf(data.X, zeros(1, size(data.X,2)), model.Sigma / model.vars(j)));
%   end
%   loss = -sum(log(P));
% catch
%   loss = inf;
% end


% d_loss.log_pis = zeros(size(state.log_pis));
% d_loss.Sigma_vec = zeros(size(state.Sigma_vec));


