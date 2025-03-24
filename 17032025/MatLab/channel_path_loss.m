function PL = channel_path_loss(PL_CO_lin, d, alpha)

% The path loss is comupted using the log-distance model
PL = PL_CO_lin*(d.^(-alpha));

end

