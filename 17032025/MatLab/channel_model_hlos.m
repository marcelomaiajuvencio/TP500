function [Hlos] = channel_model_hlos(M,N, d, kappa, N_realizacao)
% Channel Model for each device.

%% Input parameters:
% N    -> Number of Tx Antennas in the PB  (scalar)
% x, y -> (x,y) coordinates of EH devices' positions with respect to the PB
%         Each x and y can be equal-size vectors       
% kappa-> LOS factor of the Rician quasi-static fading model (scalar)

%% Output parameters:
% H    -> [N_realizacao x K x N] % Channel Matrix
% Hlos -> [N_realizacao x K x N] % LoS Components of the Channel Matrix

%% Main Code
ang = atan(d); % angular position of each EH device

% LOS Rician channel for ULA with a preventive phase shift of [0 pi 0 pi 0 ....]
for r = 1:N_realizacao
    Hlos{r} = randn(M,N) + j.*randn(M,N);
   
end

end