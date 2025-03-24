function [Hnlos] = channel_model_nlos(M,N, kappa, N_realizacao)
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
pdf1 = makedist('Rayleigh', 'b', 1);

% LOS Rician channel for ULA with a preventive phase shift of [0 pi 0 pi 0 ....]
for r = 1:N_realizacao
    %Hnlos{r} = sqrt(1/((1+kappa)))*(randn(M,N)+1i*randn(M,N)); % NLOS Rician channel
    Hnlos{r} = pdf1.random(M,N) + j.*pdf1.random(M,N);
end

end