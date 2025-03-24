% Implementação PSO - ISWCS 
% PSO - Particle Swarm Optimization
% ISWCS - International Symposium on Wireless Communication Systems
% Victoria Dala Pegorara Souto

clear all
clc
seed = 1; % Seed value
randn('seed',seed) % set the seed for random number generation - it allows reproducibility and fair comparisson of results
rng('default');
rand('seed',seed)

%% Wireless Network Parameters
N_Realizacao = 2;
N = 10; % Number of antennas at the PB
M = 20; % Number of elements at the STAR-RIS
kappa = 3.0; % Rician factor
R_min = 2; %10.^(SNR/10); 
Sigma = -80; % dBm 
Sigma_lin = (10.^(Sigma/10))./(10^3); 
Pt = 1000*(10^-3); % Total Transmit Power

% Path Loss Variables
PL_D0 = 1; % m - Reference distance
PL_C0 = -30; % dB - PL in the reference distance
PL_C0_lin = 10.^(PL_C0/10); 
dG = 70; % Distance between the BS and the STAR-RIS
dv = 2; % Vertical distance
dt = 20; % Distance horizontal between the STAR-RIS and the UE in the transmission side
dr = 50; % Distance horizontal between the STAR-RIS and the UE in the reflection side

% Channel Model
% BS -> STAR-RIS (G) - LoS channel
Alpha_BS = 2; % PL expoent
% BS -> UEr (hd) - LoS channel
Alpha_BS_UEr = 2; % PL expoent
% STAR-RIS -> UE r (hr) - Rayleight fading channel
Alpha_SUr = 2; % Expoente de PL
% STAR-RIS -> UE t (ht) - Rayleight fading channel
Alpha_SUt = 2; % Expoente de PL

% Distances
dBr = sqrt((dG - dr)^2 + dv^2);  % BS-UEr
dBt = sqrt((dG + dt)^2 + dv^2); % BS-UEt
dSt = sqrt(dt^2 + dv^2); % STAR-RIS - UEt
dSr = sqrt(dr^2 + dv^2);  % STAR-RIS - UEr

%% PSO Parameters
Nit = 10;
c1 = 0.7;
c2 = 0.7;
wmin = 0.2;
wmax = 0.9;
TamPop = 20;

%% Main Code

% Path Loss Analysis
% BS to STAR-RIS
PL_BS = channel_path_loss(PL_C0_lin, dG, Alpha_BS);

% BS to UEr
PL_BS_UEr = channel_path_loss(PL_C0_lin, dBr, Alpha_BS_UEr);

% STAR-RIS to UEr
PL_SUr = channel_path_loss(PL_C0_lin, dSr, Alpha_SUr);

% STAR-RIS to UEt
PL_SUt = channel_path_loss(PL_C0_lin, dSt, Alpha_SUt);


for real = 1:N_Realizacao
    
    fprintf('OMA - Realizacao: %d \n', real);
    
    % --> Channel Definition
 
    % BS to STAR-RIS
    hG_los = channel_model_hlos(M, N, dG, kappa, 1);
    hG_nlos = channel_model_nlos(M, N, kappa, 1);
    hG = channel_model_H(hG_los,hG_nlos,kappa,1);
    hG{1} = sqrt(PL_BS).*hG{1};
    
    % BS to UEr
    hd_los = channel_model_hlos(1, N, dBr, kappa, 1);
    hd_nlos = channel_model_nlos(1, N, kappa, 1);
    hd = channel_model_H(hd_los,hd_nlos,kappa,1);
    hd{1} = sqrt(PL_BS_UEr).*hd{1};
    
    % STAR-RIS to UEr
    hr_los = channel_model_hlos(M, 1, dSr, kappa, 1);
    hr_nlos = channel_model_nlos(M, 1, kappa, 1);
    hr = channel_model_H(hr_los,hr_nlos, kappa, 1);
    hr{1} = sqrt(PL_SUr).*hr{1};
    
    % STAR-RIS to UEt
    ht_los = channel_model_hlos(M, 1, dSt, kappa, 1);
    ht_nlos = channel_model_nlos(M, 1, kappa, 1);
    ht = channel_model_H(ht_los,ht_nlos,kappa, 1);
    ht{1} = sqrt(PL_SUt).*ht{1};

    for g = 1:Nit
        
        if (g == 1)
            
            % First Swarm
            
            % Position
            Wr_pos = exp(-1j.*2.*pi./randi(100,N,TamPop));
            Wt_pos = exp(-1j.*2.*pi./randi(100,N,TamPop));
            
            Wr_pos = (sqrt(Pt).*Wr_pos)./(sqrt(N).*abs(Wr_pos));
            Wt_pos = (sqrt(Pt).*Wt_pos)./(sqrt(N).*abs(Wt_pos));
            
            Amp_Tr_pos = rand(M,TamPop);
            Amp_Tt_pos = 1 - Amp_Tr_pos;
            
            Theta_Tr_pos = 2.*pi./randi(100,M,TamPop);
            Theta_Tt_pos = 2.*pi./randi(100,M,TamPop);
            
            Tr_pos = Amp_Tr_pos.*exp(-1j.*Theta_Tr_pos);            
            Tt_pos = Amp_Tt_pos.*exp(-1j.*Theta_Tt_pos);
            
            % Velocity
            Wr_vel = exp(-1j.*2.*pi./randi(1000,N,TamPop));
            Wt_vel = exp(-1j.*2.*pi./randi(1000,N,TamPop));
            
            Amp_Tr_vel = 0.5.*rand(M,TamPop);
            Amp_Tt_vel = 0.5.*rand(M,TamPop);
            
            Theta_Tr_vel = 2.*pi./randi(100,M,TamPop);
            Theta_Tt_vel = 2.*pi./randi(100,M,TamPop);
            
            Tr_vel = Amp_Tr_pos.*exp(-1j.*Theta_Tr_vel);
            Tt_vel = Amp_Tt_pos.*exp(-1j.*Theta_Tt_vel);
            
            power_pos(1,:) = rand(1, TamPop);
            power_pos(2,:) = 1 - power_pos(1,:);
            
            power_vel = rand(2, TamPop);
            
            % Compute the Fitness
            for pop = 1:TamPop
                SNRr(pop) = abs((hr{1}'*diag(Tr_pos(:,pop))*hG{1} + hd{1})*Wr_pos(:,pop))^2;
                SNRt(pop) = abs((ht{1}'*diag(Tt_pos(:,pop))*hG{1})*Wt_pos(:,pop))^2;
                
                Rt(pop) = 0.5*log2(1 + (Pt.*SNRt(pop)/(Sigma_lin)));
                Rr(pop) = 0.5*log2(1 + (Pt.*SNRr(pop)/(Sigma_lin)));
                
                if (Rr(pop) >= R_min && Rt(pop) >= R_min)
                    Fitness(pop) = Rr(pop) + Rt(pop);
                else
                    if (Rr(pop) < R_min)
                        if (Rt(pop) < R_min)
                            Fitness(pop) = (Rr(pop) - R_min) + (Rt(pop) - R_min);
                        else
                            Fitness(pop) = (Rr(pop) - R_min);
                        end
                    else
                        Fitness(pop) = Rt(pop) - R_min; %1000*abs(SINRt(pop) - SNR_lin);
                    end
                end
            end
            
            % Define Pbest and Gbest
            
            % Pbest
            Wr_pos_pbest = Wr_pos;
            Wt_pos_pbest = Wt_pos;
            
            power_pos_pbest = power_pos;
            
            Theta_Tr_pos_pbest = Theta_Tr_pos;
            Theta_Tt_pos_pbest = Theta_Tt_pos;
            
            Amp_Tr_pos_pbest = Amp_Tr_pos;
            Amp_Tt_pos_pbest = Amp_Tt_pos;
            
            Fitness_pbest = Fitness;
            
            % Gbest
            [Fitness_gbest pos_gbest] = max(Fitness);
            Wr_pos_gbest = Wr_pos(:,pos_gbest);
            Wt_pos_gbest = Wt_pos(:,pos_gbest);
            
            power_pos_gbest = power_pos(:, pos_gbest);
            
            Amp_Tr_pos_gbest = Amp_Tr_pos(:,pos_gbest);
            Amp_Tt_pos_gbest = Amp_Tt_pos(:,pos_gbest);
            
            Theta_Tr_pos_gbest = Theta_Tr_pos(:,pos_gbest);
            Theta_Tt_pos_gbest = Theta_Tt_pos(:,pos_gbest);
            
            % Update the Velocity and Position
            w = wmax - ((wmax - wmin)*g/Nit);
            
            for pop=1:TamPop
                Wr_vel(:,pop) = w.*Wr_vel(:,pop) + c1.*rand(1).*(Wr_pos_pbest(:,pop) - Wr_pos(:,pop)) + c2.*rand(1).*(Wr_pos_gbest - Wr_pos(:,pop));
                Wt_vel(:,pop) = w.*Wt_vel(:,pop) + c1.*rand(1).*(Wt_pos_pbest(:,pop) - Wt_pos(:,pop)) + c2.*rand(1).*(Wt_pos_gbest - Wt_pos(:,pop));

                Amp_Tr_vel(:,pop) = w.*Amp_Tr_vel(:,pop) + c1.*rand(1).*(Amp_Tr_pos_pbest(:,pop) - Amp_Tr_pos(:,pop)) + c2.*rand(1).*(Amp_Tr_pos_gbest - Amp_Tr_pos(:,pop));
                Amp_Tt_vel(:,pop) = w.*Amp_Tt_vel(:,pop) + c1.*rand(1).*(Amp_Tt_pos_pbest(:,pop) - Amp_Tt_pos(:,pop)) + c2.*rand(1).*(Amp_Tt_pos_gbest - Amp_Tt_pos(:,pop));
                
                Theta_Tr_vel(:,pop) = w.*Theta_Tr_vel(:,pop) + c1.*rand(1).*(Theta_Tr_pos_pbest(:,pop) - Theta_Tr_pos(:,pop)) + c2.*rand(1).*(Theta_Tr_pos_gbest - Theta_Tr_pos(:,pop));
                Theta_Tt_vel(:,pop) = w.*Theta_Tt_vel(:,pop) + c1.*rand(1).*(Theta_Tt_pos_pbest(:,pop) - Theta_Tt_pos(:,pop)) + c2.*rand(1).*(Theta_Tt_pos_gbest - Theta_Tt_pos(:,pop));

                power_vel(:,pop) = w.*power_vel(:,pop) + c1.*rand(1).*(power_pos_pbest(:,pop) - power_pos(:,pop)) + c2.*rand(1).*(power_pos_gbest - power_pos(:,pop));

                Theta_Tr_pos(:,pop) = Theta_Tr_pos(:,pop) + Theta_Tr_vel(:,pop);
                Theta_Tt_pos(:,pop) = Theta_Tt_pos(:,pop) + Theta_Tt_vel(:,pop);
                
                Amp_Tr_pos(:,pop) = abs(Amp_Tr_pos(:,pop) + Amp_Tr_vel(:,pop));
                Amp_Tt_pos(:,pop) = abs(Amp_Tt_pos(:,pop) + Amp_Tt_vel(:,pop));
                
                power_pos(:,pop) = power_pos(:,pop) + power_vel(:, pop);
                
                if (power_pos(1,pop) >= 1)
                    power_pos(1,pop) = rand(1);
                    power_pos(2,pop) = 1 - power_pos(1,pop);
                else
                    power_pos(2,pop) = 1 - power_pos(1,pop);
                end
                
                if (power_pos(1,pop) < 0)
                    power_pos(1,pop) = rand(1);
                    power_pos(2,pop) = 1 - power_pos(1,pop);
                else
                    power_pos(2,pop) = 1 - power_pos(1,pop);
                end
                
                for m = 1:M
                    if (Amp_Tr_pos(m,pop) < 1)
                        Amp_Tt_pos(m,pop) = 1 - Amp_Tr_pos(m,pop);
                    else
                        Amp_Tt_pos(m,pop) = rand(1);
                        Amp_Tr_pos(m,pop) = 1 - Amp_Tt_pos(m,pop);
                    end
                end
                
                Wr_pos(:,pop) = Wr_pos(:,pop) + Wr_vel(:,pop);
                Wt_pos(:,pop) = Wt_pos(:,pop) + Wt_vel(:,pop);
                
            end
            
            Wr_pos = (sqrt(Pt).*Wr_pos)./(sqrt(N).*abs(Wr_pos));
            Wt_pos = (sqrt(Pt).*Wt_pos)./(sqrt(N).*abs(Wt_pos));
            
            Tr_pos = Amp_Tr_pos.*exp(j.*Theta_Tr_pos);            
            Tt_pos = Amp_Tt_pos.*exp(j.*Theta_Tt_pos);
            
            % Compute the Fitness
            for pop = 1:TamPop
                SNRr(pop) = abs((hr{1}'*diag(Tr_pos(:,pop))*hG{1} + hd{1})*Wr_pos(:,pop))^2;
                SNRt(pop) = abs((ht{1}'*diag(Tt_pos(:,pop))*hG{1})*Wt_pos(:,pop))^2;
                
                Rt(pop) = 0.5*log2(1 + (Pt.*SNRt(pop)/(Sigma_lin)));
                Rr(pop) = 0.5*log2(1 + (Pt.*SNRr(pop)/(Sigma_lin)));
                
                if (Rr(pop) >= R_min && Rt(pop) >= R_min)
                    Fitness(pop) = Rr(pop) + Rt(pop);
                else
                    if (Rr(pop) < R_min)
                        if (Rt(pop) < R_min)
                            Fitness(pop) = (Rr(pop) - R_min) + (Rt(pop) - R_min);
                        else
                            Fitness(pop) = (Rr(pop) - R_min);
                        end
                    else
                        Fitness(pop) = Rt(pop) - R_min; %1000*abs(SINRt(pop) - SNR_lin);
                    end
                end
            end
            
            Best(real,g) = Fitness_gbest;
            
            
            
        else % (g != 1)
            
            % Define Pbest and Gbest
            
            % Pbest
            for pop = 1:TamPop
                if (Fitness_pbest(pop) < Fitness(pop))
                    Wr_pos_pbest(:,pop) = Wr_pos(:,pop);
                    Wt_pos_pbest(:,pop) = Wt_pos(:,pop);
                    
                    power_pos_pbest(:,pop) = power_pos(:,pop);
                    
                    Theta_Tr_pos_pbest(:,pop) = Theta_Tr_pos(:,pop);
                    Theta_Tt_pos_pbest(:,pop) = Theta_Tt_pos(:,pop);

                    Amp_Tr_pos_pbest(:,pop) = Amp_Tr_pos(:,pop);
                    Amp_Tt_pos_pbest(:,pop) = Amp_Tt_pos(:,pop);

                    Fitness_pbest(pop) = Fitness(pop);
                end
            end
            
            % Gbest
            if (Fitness_gbest < max(Fitness))
                [Fitness_gbest pos_gbest] = max(Fitness);
                Wr_pos_gbest = Wr_pos(:,pos_gbest);
                Wt_pos_gbest = Wt_pos(:,pos_gbest);
                
                power_pos_gbest = power_pos(:, pos_gbest);
                
                Amp_Tr_pos_gbest = Amp_Tr_pos(:,pos_gbest);
                Amp_Tt_pos_gbest = Amp_Tt_pos(:,pos_gbest);

                Theta_Tr_pos_gbest = Theta_Tr_pos(:,pos_gbest);
                Theta_Tt_pos_gbest = Theta_Tt_pos(:,pos_gbest);
            end
            
            % Update the Velocity and Position
            w = wmax - ((wmax - wmin)*g/Nit);
            
            for pop=1:TamPop
                
                Wr_vel(:,pop) = w.*Wr_vel(:,pop) + c1.*rand(1).*(Wr_pos_pbest(:,pop) - Wr_pos(:,pop)) + c2.*rand(1).*(Wr_pos_gbest - Wr_pos(:,pop));
                Wt_vel(:,pop) = w.*Wt_vel(:,pop) + c1.*rand(1).*(Wt_pos_pbest(:,pop) - Wt_pos(:,pop)) + c2.*rand(1).*(Wt_pos_gbest - Wt_pos(:,pop));
                
                Amp_Tr_vel(:,pop) = w.*Amp_Tr_vel(:,pop) + c1.*rand(1).*(Amp_Tr_pos_pbest(:,pop) - Amp_Tr_pos(:,pop)) + c2.*rand(1).*(Amp_Tr_pos_gbest - Amp_Tr_pos(:,pop));
                Amp_Tt_vel(:,pop) = w.*Amp_Tt_vel(:,pop) + c1.*rand(1).*(Amp_Tt_pos_pbest(:,pop) - Amp_Tt_pos(:,pop)) + c2.*rand(1).*(Amp_Tt_pos_gbest - Amp_Tt_pos(:,pop));
                
                Theta_Tr_vel(:,pop) = w.*Theta_Tr_vel(:,pop) + c1.*rand(1).*(Theta_Tr_pos_pbest(:,pop) - Theta_Tr_pos(:,pop)) + c2.*rand(1).*(Theta_Tr_pos_gbest - Theta_Tr_pos(:,pop));
                Theta_Tt_vel(:,pop) = w.*Theta_Tt_vel(:,pop) + c1.*rand(1).*(Theta_Tt_pos_pbest(:,pop) - Theta_Tt_pos(:,pop)) + c2.*rand(1).*(Theta_Tt_pos_gbest - Theta_Tt_pos(:,pop));

                power_vel(:,pop) = w.*power_vel(:,pop) + c1.*rand(1).*(power_pos_pbest(:,pop) - power_pos(:,pop)) + c2.*rand(1).*(power_pos_gbest - power_pos(:,pop));

                Theta_Tr_pos(:,pop) = Theta_Tr_pos(:,pop) + Theta_Tr_vel(:,pop);
                Theta_Tt_pos(:,pop) = Theta_Tt_pos(:,pop) + Theta_Tt_vel(:,pop);
                
                Amp_Tr_pos(:,pop) = abs(Amp_Tr_pos(:,pop) + Amp_Tr_vel(:,pop));
                Amp_Tt_pos(:,pop) = abs(Amp_Tt_pos(:,pop) + Amp_Tt_vel(:,pop));
                
                power_pos(:,pop) = power_pos(:,pop) + power_vel(:, pop);
                
                if (power_pos(1,pop) >= 1)
                    power_pos(1,pop) = rand(1);
                    power_pos(2,pop) = 1 - power_pos(1,pop);
                else
                    power_pos(2,pop) = 1 - power_pos(1,pop);
                end
                
                if (power_pos(1,pop) < 0)
                    power_pos(1,pop) = rand(1);
                    power_pos(2,pop) = 1 - power_pos(1,pop);
                else
                    power_pos(2,pop) = 1 - power_pos(1,pop);
                end
                                
                for m = 1:M
                    if (Amp_Tr_pos(m,pop) < 1)
                        Amp_Tt_pos(m,pop) = 1 - Amp_Tr_pos(m,pop);
                    else
                        Amp_Tt_pos(m,pop) = rand(1);
                        Amp_Tr_pos(m,pop) = 1 - Amp_Tt_pos(m,pop);
                    end
                end
                
                
                Wr_pos(:,pop) = Wr_pos(:,pop) + Wr_vel(:,pop);
                Wt_pos(:,pop) = Wt_pos(:,pop) + Wt_vel(:,pop);
               
            end
            
            Wr_pos = (sqrt(Pt).*Wr_pos)./(sqrt(N).*abs(Wr_pos));
            Wt_pos = (sqrt(Pt).*Wt_pos)./(sqrt(N).*abs(Wt_pos));
            
            
            Tr_pos = Amp_Tr_pos.*exp(j.*Theta_Tr_pos);            
            Tt_pos = Amp_Tt_pos.*exp(j.*Theta_Tt_pos);
            
            % Compute the Fitness
            for pop = 1:TamPop
                SNRr(pop) = abs((hr{1}'*diag(Tr_pos(:,pop))*hG{1} + hd{1})*Wr_pos(:,pop))^2;
                SNRt(pop) = abs((ht{1}'*diag(Tt_pos(:,pop))*hG{1})*Wt_pos(:,pop))^2;
                
                Rt(pop) = 0.5*log2(1 + (Pt.*SNRt(pop)/(Sigma_lin)));
                Rr(pop) = 0.5*log2(1 + (Pt.*SNRr(pop)/(Sigma_lin)));
                
                if (Rr(pop) >= R_min && Rt(pop) >= R_min)
                    Fitness(pop) = Rr(pop) + Rt(pop);
                else
                    if (Rr(pop) < R_min)
                        if (Rt(pop) < R_min)
                            Fitness(pop) = (Rr(pop) - R_min) + (Rt(pop) - R_min);
                        else
                            Fitness(pop) = (Rr(pop) - R_min);
                        end
                    else
                        Fitness(pop) = Rt(pop) - R_min; %1000*abs(SINRt(pop) - SNR_lin);
                    end
                end
            end
            
        end
        
        if (Fitness_gbest ~= 0)
            Best(real,g) = Fitness_gbest;
        else
            Best(real,g) = 0;
        end
    end
    
end
%Fitness_gbest
%Sol = (mean(Best));
%mean(Best)
Sol = mean(Best(:));
Sol
%disp(Best(2,2));
%Best


