%% ver5 
% ITNTN Simulation
% Multi-User Integrated Terrestrial and Non-Terrestrial Network

%% Clears workspace, command window and closes all figures
clear; clc; close all;

%% Simulation Parameters
% Altitude of LEO satellite (m)
H = 450 * 10^3;
% Coordinate of LEO satellite (m)
s_x = 0; s_y = 0; s_z = H;      

% Coverage radius of BS (m)
D = 2000;                       
% Coordinates of BS (m)
b_x = 0; b_y = 0; b_z = 0;

% Carrier frequency (C-band) (Hz)
f_Ter_ac = 6 * 10^9;             

% Channel Bandwidths (Hz)
W_Ter = 5 * 10^6; W_Sat = 30 * 10^6;                  

% Transmission Powers (dBm → W)
P_Ter_tx = 10^((42-30)/10); % 42 dBm → ~15.85 W
P_Sat_tx = 10^((60-30)/10); % 60 dBm → 1000 W    

% Noise Figures (dB)
N_Ter = 0.8; N_Sat = 0.8;  
% Noise Spectral Density (W/Hz) 
% N0 = kB * T0 * NF where kB = 1.38 * 10^-23 J/K and T0 = 290K
N0_Ter = 1.38e-23 * 290 * 10^(N_Ter/10);
N0_Sat = 1.38e-23 * 290 * 10^(N_Sat/10);

% Path loss exponent
eta = 2.2;

% Base prices
alpha_Ter = 0.5; alpha_Sat = 0.85;

% Load sensitivity factors
beta_Ter = 1.5; beta_Sat = 1.2;

% Available capacity sensitivity factors
gamma_Ter = 0.05; gamma_Sat = 0.03;

% Energy to-cost conversion factor
theta = 1.5; 

% Total cost weights
% C_T = xi_a * C_A + xi_b * C_E
% xi_a + xi_b = 1
% Equal weighting between network access and energy consumption costs
xi_a = 0.5; xi_b = 0.5;

%% Network Capacity
% Maximum capacity (bbu)
C_Ter = 45; C_Sat = 50;

% Admission thresholds for new calls (bbu)
t_m_Ter = 35; t_m_Sat = 30;

% Reserves capacity for handoff calls (pioritises handoff calls)
% Reserved capacity = Maximum capacity - Admission threshold
 
fprintf('Network Capacity:\n');
fprintf('  Terrestrial: C = %d bbu, t_m = %d bbu\n', C_Ter, t_m_Ter);
fprintf('  Satellite:   C = %d bbu, t_m = %d bbu\n', C_Sat, t_m_Sat);

%% Initial Network State
% Initial network loads (bbu)
L_Ter = 5; L_Sat = 5;

% Available capacity (%)
C_Ter_avail = 1 - (L_Ter / C_Ter); C_Sat_avail = 1 - (L_Sat / C_Sat);

fprintf('\nInitial Network State:\n');
fprintf(['  Terrestrial: Network Load = %.2f%%, ' ...
         'Available Capacity = %.2f%%\n'], ...
    (L_Ter/C_Ter) * 100, C_Ter_avail * 100);
fprintf(['  Satellite:   Network Load = %.2f%%,' ...
         'Available Capacity = %.2f%%\n'], ...
    (L_Sat/C_Sat) * 100, C_Sat_avail * 100);

%% Multi-User Configuration
% Number of users
U = 20;
% Number of position updates
K = 16;                         

% File size (bits)
file_size = 1e6;
% Allocated bandwidth per call (bbu)
b_user = 2;  

% Assigns same file size and allocated bandwidth to all users
O_k = file_size * ones(1, U); b_k = b_user * ones(1, U);  

%% User Mobility Paths
user_paths = cell(U, 1);
for u = 1:U
    if mod(u, 2) == 1
        % Odd users: Start within BS coverage, move outward
        ux_path = linspace(0, 3000, K);    
        uy_path = (u - 1) * 50 + 200;
    else
        % Even users: Start outside BS coverage, move inward
        ux_path = linspace(3000, 0, K);    
        uy_path = (u - 1) * 50 + 250;
    end
    user_paths{u} = [ux_path; uy_path * ones(1,K); zeros(1,K)];
end

%% QoS Parameters
% New call arrival rates
lambda_m_Ter = 2.0; lambda_m_Sat = 1.8;
% Handoff call arrival rates
lambda_h_Ter = 0.8; lambda_h_Sat = 0.6;  

% New call departure rates
mu_m_Ter = 1; mu_m_Sat = 1;
% Handoff call departure rates
mu_h_Ter = 0.5; mu_h_Sat = 0.5;          

%% Result Arrays
% Stores network selection 
% 0 = blocked, 1 = Ter, 2 = Sat, -1 = dropped
all_decisions = zeros(U, K);

% Stores total cost for selected network
all_costs = zeros(U, K);
% Stores per-network total costs
all_costs_Ter = zeros(U, K); all_costs_Sat = zeros(U, K); 

% Stores the coverage status 
% 1 = Within BS coverage, 0 = outside Ter coverage)
all_coverage = zeros(U, K);

% Stores per-network admission status
% 'feasible', 'blocked', 'dropped'
all_status_Ter = cell(U, K); all_status_Sat = cell(U, K);

% Initialises number of handoff calls
handoff_calls = 0;

% Initialises number of dropped calls
dropped_call = 0;

% Initialises number of blocked calls
blocked_calls = 0;

%% Multi-User Simulation Loop
fprintf('\nSimulation: \n');

% Outer loop: Iterates over users
for u = 1:U
    % Gets current user's file size (bits)   
    O = O_k(u); 
    % Gets current user's allocated Bandwidth (bbu)
    b_user = b_k(u);
    % Gets currrent user's mobility path
    path = user_paths{u};
    
    % Initialises previous network selection
    X_prev = [];

    % Initialises flag to check if call was blocked or dropped
    isCallBlockedOrDropped = false;
    
    % Inner loop: Iterates over position updates
    for k = 1:K
        % Get current user position
        u_x = path(1,k); u_y = path(2,k); u_z = path(3,k);
        
        %% Terrestrial Access
        % Computes communication distance between user and BS
        l_Ter_ac = sqrt((u_x - b_x)^2 + (u_y - b_y)^2 + (u_z - b_z)^2);
        
        % Checks coverage status
        % Return 1 if user within BS coverage, else 0
        all_coverage(u,k) = (l_Ter_ac <= D);
        
        % Computes path loss 
        PL_dB = 52.44 + 20 * log10(l_Ter_ac / 10^3) + 20 * log10(f_Ter_ac / 10^6);
        
        % Generates complex Rayleigh fading coefficient (g ~ CN(0,1))
        g = (randn(1) + 1i * randn(1)) / sqrt(2);

        % Computes Channel gain
        G_Ter = abs(g)^2 / 10^(PL_dB/10);
        
        % Computes achievable transmission rate (Shannon Capacity Theorem)
        R_Ter_tx = W_Ter * log2(1 + (P_Ter_tx * G_Ter)/(W_Ter * N0_Ter));
        
        % Computes transmission time
        T_Ter_tx = O / R_Ter_tx;
        
        %% Satellite Access
        % Compuutes communication distance between user and LEO satellite
        l_Sat_ac = sqrt((u_x - s_x)^2 + (u_y - s_y)^2 + (u_z - s_z)^2);
        
        % Computes channel gain
        G_Sat = (l_Sat_ac / 10^3)^(-eta);
        
        % Computes achievable transmission rate (Shannon Capacity Theorem)
        R_Sat_tx = W_Sat * log2(1 + (P_Sat_tx * G_Sat)/(W_Sat * N0_Sat));
        
        % Computes transmission time
        T_Sat_tx = O / R_Sat_tx;
        
        % Initialises total costs
        C_T_Ter = 0; C_T_Sat = 0;
        
        %% At k = 1, New calls
        if k == 1
            % Checks if networks can admit new calls
            % Returns admission status and updated network load and available capacity
            [status_Ter, L_Ter, C_Ter_avail] = ...
                check_admission(L_Ter, C_Ter_avail, b_user, C_Ter, t_m_Ter, 'new');
            
            [status_Sat, L_Sat, C_Sat_avail] = ...
                check_admission(L_Sat, C_Sat_avail, b_user, C_Sat, t_m_Sat, 'new');
            
            % Compute costs if feasible
            if strcmp(status_Ter, 'feasible')
                % Computes unit price (based on network load and available capacity)
                p_Ter = alpha_Ter + beta_Ter * (L_Ter / C_Ter) + gamma_Ter / C_Ter_avail;
                
                % Computes network access cost
                C_A_Ter = p_Ter * b_user;
                
                % Computes energy consumption
                E_Ter = P_Ter_tx * T_Ter_tx;
                
                % Computes energy cost
                C_E_Ter = theta * E_Ter;
                
                % Computes total cost
                C_T_Ter = xi_a * C_A_Ter + xi_b * C_E_Ter;
            end
            
            if strcmp(status_Sat, 'feasible')
                % Computes unit price (based on network load and available capacity)
                p_Sat = alpha_Sat + beta_Sat * (L_Sat/C_Sat) + gamma_Sat / C_Sat_avail;
                
                % Computes network access cost
                C_A_Sat = p_Sat * b_user;
                
                % Computes energy consumption
                E_Sat = P_Sat_tx * T_Sat_tx;
                % Computes Eenergy consumption cost
                C_E_Sat = theta * E_Sat;
                
                % Computes total cost
                C_T_Sat = xi_a * C_A_Sat + xi_b * C_E_Sat;
            end
            
            % Store per-network total costs
            all_costs_Ter(u,k) = C_T_Ter; all_costs_Sat(u,k) = C_T_Sat;
            
            %% Network Selection
            % If the user is outside BS coverage
            if (all_coverage(u,k) == 0)
                % If the satellite network is feasible
                if (strcmp(status_Sat, 'feasible'))
                    X = [0; 1]; selected_network = 2;
                % If the satellite network is not feasible → blocks the call
                else
                    X = [0; 0]; selected_network = 0;
                    
                    % Increments number of blocked calls
                    blocked_calls = blocked_calls + 1;

                    % Sets the flag to true
                    isCallBlockedOrDropped = true;
                end
            end    
            
            % If the user is within BS coverage
            if (all_coverage(u, k) == 1)
                % Case 1: if both networks are feasible
                if (strcmp(status_Ter, 'feasible') && ...
                    strcmp(status_Sat, 'feasible'))
                    % Compares total costs
                    if C_T_Ter < C_T_Sat
                        X = [1; 0]; selected_network = 1; 
                    else
                        X = [0; 1]; selected_network = 2; 
                    end
                % Case 2: If only terrestrial network is feasible
                elseif (strcmp(status_Ter, 'feasible'))
                    X = [1; 0]; selected_network = 1; 
                
                % Case 3: If only satellite network is feasible:
                elseif (strcmp(status_Sat, 'feasible'))
                    X = [0; 1]; selected_network = 2; 
                
                % Case 4: If both networks are not feasible → block the call
                else    
                    X = [0; 0]; selected_network = 0;
                    
                    % Increments number of blocked calls
                    blocked_calls = blocked_calls + 1;

                    % Sets the flag to true
                    isCallBlockedOrDropped = true;
                end
            end

            % Stores the selected network
            all_decisions(u,k) = selected_network;

            % Stores the selected network total cost
            all_costs(u,k) = X' * [C_T_Ter; C_T_Sat];

            % Stores the admission status
            all_status_Ter{u,k} = status_Ter; 
            all_status_Sat{u,k} = status_Sat;

            % Updates previous network selection
            X_prev = X;
            
        %% At k = 2..K, Handoff calls
        else
            % If the call was blocked or dropped, skip remaining position updates
            % Carries forward previous network selection decision (B or D)
            if isCallBlockedOrDropped
                all_decisions(u, k) = all_decisions(u, k-1);
                all_costs(u,k) = 0;
                continue;
            end
            
            % Terrestrial:
            % If the user is within terrestrial coverage
            if (all_coverage(u,k) == 1)
                % If the terrestrial network is feasible
                if (strcmp(all_status_Ter{u,1}, 'feasible'))
                    % Computes unit price (based on network load and available capacity)
                    p_Ter = alpha_Ter + beta_Ter * (L_Ter / C_Ter) + gamma_Ter / C_Ter_avail;
                    
                    % Computes network access cost
                    C_A_Ter = p_Ter * b_user;
                    
                    % Computes energy consumption 
                    E_Ter = P_Ter_tx * T_Ter_tx;
                    % Computes energy consumption cost
                    C_E_Ter = theta * E_Ter;
                    
                    % Computes total cost 
                    C_T_Ter = xi_a * C_A_Ter + xi_b * C_E_Ter;
                end     
            end
            
            % Satellite:
            % if the satellite is feasible
            if (strcmp(all_status_Sat{u,1}, 'feasible'))
                % Computes unit price (based on network load and available capacity)
                p_Sat = alpha_Sat + beta_Sat * (L_Sat / C_Sat) + gamma_Sat / C_Sat_avail;
                % Computes network access cost
                C_A_Sat = p_Sat * b_user;
                
                % Computes energy consumption
                E_Sat = P_Sat_tx * T_Sat_tx;
                % Computes energy consumption cost
                C_E_Sat = theta * E_Sat;
                
                % Computes total cost 
                C_T_Sat = xi_a * C_A_Sat + xi_b * C_E_Sat;
            end
            
            % Stores per-network total costs
            all_costs_Ter(u,k) = C_T_Ter; all_costs_Sat(u,k) = C_T_Sat;
            
            %% Network Selection
            % Initialise current network selection
            X = [];
            
            % If the user is within terrestrial coverage
            if (all_coverage(u,k) == 1)
                % Case 1: If both networks are feasible
                if (strcmp(all_status_Ter{u,1}, 'feasible') && ...
                    strcmp(all_status_Sat{u,1}, 'feasible'))
                    % Compare total costs
                    if (C_T_Ter < C_T_Sat)
                        X = [1; 0];
                    else
                        X = [0; 1];
                    end
                
                % Case 2: If only terrestrial network is feasible
                elseif (strcmp(all_status_Ter{u,1}, 'feasible'))
                    X = [1; 0];
                
                % Case 3: If only satellite network is feasible
                elseif (strcmp(all_status_Sat{u,1}, 'feasible'))
                   X = [0; 1];
                end
            
            % if the user is outside BS coverage
            else
                % if the satellite network is feasible
                if (strcmp(all_status_Sat{u,1}, 'feasible'))
                    X = [0; 1];
                end
            end
            
            %% Handoff Execution
            % Checks if the current network selection (X) is the same as
            % previous network selection (X_prev)
            if isequal(X_prev, X)
                % Yes: Stay on the same network (No handoff needed)
                selected_network = all_decisions(u, k-1);
            else
                % No: Attempt handoff
                % Handoff to Terrestrial
                if isequal(X, [1; 0])
                    % Checks admission for handoff call
                    [status_Ter, L_Ter_temp, C_Ter_avail_temp] = ...
                        check_admission(L_Ter, C_Ter_avail, b_user, ...
                        C_Ter, t_m_Ter, 'handoff');
                    
                    % If the terrestrial network is feasible → Handoff successful
                    if (strcmp(status_Ter, 'feasible'))
                        % Updates network state
                        L_Ter = L_Ter_temp; 
                        C_Ter_avail = C_Ter_avail_temp;
                       
                        % Releases resources from previous network
                        L_Sat = L_Sat - b_user; 
                        C_Sat_avail = 1 - (L_Sat / C_Sat);
                        
                        selected_network = 1;
                        
                        % Increments the number of handoff calls
                        handoff_calls = handoff_calls + 1;
                    
                    % If the terrestrial network is not feasible → Handoff failed
                    else
                        X = [0; 0]; selected_network = -1; 
                        
                        % Increments number of dropped calls
                        dropped_call = dropped_call + 1;
                        
                        % Sets the flag to true
                        isCallBlockedOrDropped = true;
                        
                        % Releases resources from previous network
                        if isequal(X_prev, [1; 0])
                            L_Ter = L_Ter - b_user; 
                            C_Ter_avail = 1 - (L_Ter / C_Ter);
                        elseif isequal(X_prev, [0; 1])
                            L_Sat = L_Sat - b_user; 
                            C_Sat_avail = 1 - (L_Sat / C_Sat);
                        end
                    end
                
                % Handoff to Satellite
                elseif isequal(X, [0; 1])
                    % Checks admission for handoff call
                    [status_Sat, L_Sat_temp, C_Sat_avail_temp] = ...
                        check_admission(L_Sat, C_Sat_avail, b_user, ...
                        C_Sat, t_m_Sat, 'handoff');
                    
                    % If the satellite network is feasible → Handoff successful
                    if (strcmp(status_Sat, 'feasible'))
                        % Updates network state
                        L_Sat = L_Sat_temp; 
                        C_Sat_avail = C_Sat_avail_temp;
                        
                        % Releases resources from previous network
                        L_Ter = L_Ter - b_user; 
                        C_Ter_avail = 1 - (L_Ter / C_Ter);
                                                
                        selected_network = 2;
                        
                        % Increments number of handoff calls
                        handoff_calls = handoff_calls + 1;
                    else
                        % If the terrestrial network is not feasible → Handoff failed
                        % 
                        X = [0; 0]; selected_network = -1; 

                        % Increments the number of dropped calls
                        dropped_call = dropped_call + 1;
                        
                        % Sets the flag to true
                        isCallBlockedOrDropped = true;

                        % Releasess resources from previous network
                        if (isequal(X_prev, [1; 0]))
                            L_Ter = L_Ter - b_user; 
                            C_Ter_avail = 1 - (L_Ter / C_Ter);
                        elseif (isequal(X_prev, [0; 1]))
                            L_Sat = L_Sat - b_user; 
                            C_Sat_avail = 1 - (L_Sat / C_Sat);
                        end
                    end
                
                % If both networks are not feasible → Handoff failed    
                else
                    X = [0; 0]; selected_network = -1; 
                    
                    % Increments the number of dropped calls
                    dropped_call = dropped_call + 1;
                    
                    isCallBlockedOrDropped = true;
                    % Releases resources from previous network
                    if (isequal(X_prev, [1; 0]))
                        L_Ter = L_Ter - b_user; 
                        C_Ter_avail = 1 - (L_Ter / C_Ter);
                    elseif (isequal(X_prev, [0; 1]))
                        L_Sat = L_Sat - b_user; 
                        C_Sat_avail = 1 - (L_Sat / C_Sat);
                    end
                end
            end
            
            % Stores the selected network
            all_decisions(u,k) = selected_network;
            
            % Store the selected network total cost
            all_costs(u,k) = X' * [C_T_Ter; C_T_Sat];
            
            % Updates previous network selection
            X_prev = X;
        end
    end
end

%% Dynamic Call Blocking and Dropping Probabilities (Subplots)
% Compute probabilities per position update
P_B_over_k = zeros(1,K); P_D_over_k = zeros(1,K);

for k = 1:K
    if k == 1
        % New calls
        total_new_calls = U;
        blocked_calls_k = nnz(all_decisions(:,k) == 0);
        P_B_over_k(k) = blocked_calls_k / total_new_calls;
        P_D_over_k(k) = 0;
    else
        % Handoff calls
        total_handoff_k = 0;
        dropped_calls_k = 0;
        for u = 1:U
            prev_net = all_decisions(u,k-1);
            curr_net = all_decisions(u,k);
            if prev_net > 0
                total_handoff_k = total_handoff_k + 1;
                if curr_net == -1
                    dropped_calls_k = dropped_calls_k + 1;
                end
            end
        end
        P_B_over_k(k) = 0;
        P_D_over_k(k) = dropped_calls_k / (total_handoff_k + eps);
    end
end

%% Network Operator Revenue Per Position Update and Cumulative
% Computes the operator revenue per position update
R_over_time = sum(all_costs, 1);  

% Computes the cumulative operator revenue
R_cumulative = cumsum(R_over_time);

% %% QoS Performance Metrics: Call Blocking and Dropping Probability
% 
% % Markov Model
% % Arrival rate for new calls
% lambda_m_vals = linspace(0.5, 15, 25);
% 
% % Arrival rate for handoff calls
% lambda_h_vals = 0.5 * lambda_m_vals;                  
% 
% % Initialise probabilities
% P_B_Ter = zeros(size(lambda_m_vals));
% P_D_Ter = zeros(size(lambda_m_vals));
% P_B_Sat = zeros(size(lambda_m_vals));
% P_D_Sat = zeros(size(lambda_m_vals));
% 
% for i = 1:length(lambda_m_vals)
%     % Arrival rate for new calls
%     lambda_m_Ter = lambda_m_vals(i); lambda_m_Sat = lambda_m_vals(i);
% 
%     % Arrival rate for handoff calls
%     lambda_h_Ter = lambda_h_vals(i); lambda_h_Sat = lambda_h_vals(i);
% 
%     % Network load factor for new calls
%     rho_m_Ter = lambda_m_Ter / mu_m_Ter; 
%     rho_m_Sat = lambda_m_Sat / mu_m_Sat; 
% 
%     % Network load factor for handoff calls
%     rho_h_Ter = lambda_h_Ter / mu_h_Ter;
%     rho_h_Sat = lambda_h_Sat / mu_h_Sat;
% 
%     % Initialise state probabilities
%     ST_Ter = 0; SB_Ter = 0; SD_Ter = 0;
%     ST_Sat = 0; SB_Sat = 0; SD_Sat = 0;
% 
%     % Terrestrial Network
%     for m = 0:t_m_Ter
%         for h = 0:C_Ter
%             if (b_user * m <= t_m_Ter) && (b_user * (m + h) <= C_Ter)
% 
%                 P = (rho_m_Ter^(m) * rho_h_Ter^(h)) / (factorial(m) * ...
%                      factorial(h));
% 
%                 ST_Ter = ST_Ter + P;
% 
%                 % Blocking state: new call rejected
%                 if ((b_user * m  > t_m_Ter) || ...
%                     (b_user + b_user * (m + h) > C_Ter))
%                     SB_Ter = SB_Ter + P;
%                 end
% 
%                 % Dropping state: handoff rejected
%                 if (b_user + b_user * (m + h) > C_Ter)
%                     SD_Ter = SD_Ter + P;
%                 end
%             end
%         end
%     end
% 
%     % Satellite Network 
%     for m = 0:t_m_Sat
%         for h = 0:C_Sat
%             if (b_user * m <= t_m_Sat) && (b_user * (m + h) <= C_Sat)
%                 P = (rho_m_Sat^(m)  * rho_h_Sat^(h)) / (factorial(m) * ...
%                      factorial(h));
%                 ST_Sat = ST_Sat + P;
%                 if ((b_user*m > t_m_Sat) || ...
%                     (b_user +  b_user * (m + h) > C_Sat))
%                     SB_Sat = SB_Sat + P;
%                 end
%                 if (b_user +  b_user * (m + h) > C_Sat)
%                     SD_Sat = SD_Sat + P;
%                 end
%             end
%         end
%     end
% 
%     % Compute Probabilities
%     P_B_Ter(i) = SB_Ter / ST_Ter;  P_D_Ter(i) = SD_Ter / ST_Ter;
%     P_B_Sat(i) = SB_Sat / ST_Sat; P_D_Sat(i) = SD_Sat / ST_Sat;
% end

%% Network Operator Revenue
% Total cost per user
total_user_costs = sum(all_costs, 2);
% Network operator revenue
R_total = sum(total_user_costs);            

% Count network selections
ter_selections = nnz(all_decisions(:) == 1);
sat_selections = nnz(all_decisions(:) == 2);

%% Network Selection History Per User
for u = 1:U
    
    % Count network usage
    user_ter = nnz(all_decisions(u,:) == 1);
    user_sat = nnz(all_decisions(u,:) == 2);
    user_blocked = nnz(all_decisions(u,:) == 0);
    user_dropped = nnz(all_decisions(u,:) == -1);
   
    % Count number of handoffs for current user
    user_handoffs = 0;
    for k = 2:K
        if (all_decisions(u, k) ~= all_decisions(u, k - 1) && ...
            all_decisions(u, k) > 0 && all_decisions(u, k - 1) > 0)
            user_handoffs = user_handoffs + 1;
        end
    end
end

%% Network Operator Dashboard

% Create the main dashboard figure
dashboardFig = uifigure('Name', 'Network Operator Dashboard', ...
                        'Position', [100, 100, 1400, 900]);

% Create tab group
tabGroup = uitabgroup(dashboardFig, 'Position', [10, 10, 1380, 880]);

%% Home Tab - 3D Network Setup Visualisation

homeTab = uitab(tabGroup, 'Title', 'Home');

% Convert to km for visualisation
km_scale = 1000;

%% Base Station
% Create panel for BS 3D view
bs3DPanel = uipanel(homeTab, 'Title', 'Base Station 3D View', ...
                    'Position', [10, 440, 680, 390]);
ax_bs_3d = uiaxes(bs3DPanel, 'Position', [50, 30, 600, 320]);

% Plot BS 3D view
hold(ax_bs_3d, 'on'); grid(ax_bs_3d, 'on');

% Create circular disk for coverage area at ground level
theta_circle = linspace(0, 2*pi, 100);
[X_disk, Y_disk] = meshgrid(linspace(-D/km_scale, D/km_scale, 50), ...
                             linspace(-D/km_scale, D/km_scale, 50));
Z_disk = zeros(size(X_disk));

R_disk = sqrt(X_disk.^2 + Y_disk.^2);
Z_disk(R_disk > D/km_scale) = NaN;

% Plot coverage disk
surf(ax_bs_3d, X_disk, Y_disk, Z_disk, 'FaceColor', [0.7, 1, 0.7], ...
     'FaceAlpha', 0.5, 'EdgeColor', [0, 0.5, 0], 'EdgeAlpha', 0.3, ...
     'DisplayName', 'Coverage Area');

% Plot base station
plot3(ax_bs_3d, b_x/km_scale, b_y/km_scale, b_z/km_scale + 0.05, ...
      'r^', 'MarkerSize', 20, 'MarkerFaceColor', 'r', 'DisplayName', ...
      'Base Station');

xlabel(ax_bs_3d, 'X (km)'); ylabel(ax_bs_3d, 'Y (km)'); zlabel(ax_bs_3d, ...
      'Z (km)');
title(ax_bs_3d, sprintf('Base Station 3D View'));
legend(ax_bs_3d, 'Location', 'northeast');
view(ax_bs_3d, -37.5, 30);
axis(ax_bs_3d, 'equal');
xlim(ax_bs_3d, [-3, 3]); ylim(ax_bs_3d, [-3, 3]); zlim(ax_bs_3d, [0, 0.5]);
hold(ax_bs_3d, 'off');

% Create text panel for BS information
bsInfoPanel = uipanel(homeTab, 'Title', 'Base Station Information', ...
                      'Position', [10, 10, 680, 420]);
bsInfoText = uitextarea(bsInfoPanel, 'Position', [20, 20, 640, 380]);
bsInfoText.Editable = 'off';
bsInfoText.FontSize = 11;

% Build BS information string
bsInfoString = sprintf(['Location:\n' ...
    '  ▬ Position: (%.1f, %.1f, %.1f) km\n' ...
    '  ▬ Coverage Radius: %.1f km\n\n' ...
    ...
    'Technical Specifications:\n' ...
    '  ▬ Carrier Frequency: %.1f GHz\n' ...
    '  ▬ Channel Bandwidth: %.1f MHz\n' ...
    '  ▬ Transmission Power: %.2f W (%.1f dBm)\n' ...
    '  ▬ Noise Figure: %.1f dB\n' ...
    ...
    '\nNetwork Capacity:\n' ...
    '  ▬ Maximum Network Capacity: %d bbu\n' ...
    '  ▬ Admission Threshold (New Calls): %d bbu\n' ...
    '  ▬ Reserved Capacity (Handoff Calls): %d bbu\n\n' ...
    ...
    'Pricing Parameters:\n' ...
    '  ▬ Base Price (α): %.2f\n' ...
    '  ▬ Network Load Sensitivity (β): %.2f\n' ...
    '  ▬ Available Capacity Sensitivity (γ): %.3f'], ...
    b_x / km_scale, b_y / km_scale, b_z / km_scale, D / km_scale, ...
    f_Ter_ac / 10^9, W_Ter / 10^6, P_Ter_tx, 10 * log10(P_Ter_tx * 1000), ...
    N_Ter, ...
    C_Ter, t_m_Ter, C_Ter - t_m_Ter, ...
    alpha_Ter, beta_Ter, gamma_Ter);

bsInfoText.Value = bsInfoString;

%% LEO Satellite

% Create panel for LEO Satellite 3D view
sat3DPanel = uipanel(homeTab, 'Title', 'LEO Satellite 3D View', ...
                     'Position', [700, 440, 680, 390]);
ax_sat_3d = uiaxes(sat3DPanel, 'Position', [50, 30, 600, 320]);

% Plot LEO satellite 3D view
hold(ax_sat_3d, 'on');
grid(ax_sat_3d, 'on');

% Plot ground plane (reference) - smaller scale
ground_size = 100; % km
[X_ground, Y_ground] = meshgrid(linspace(-ground_size, ground_size, 20), ...
                                 linspace(-ground_size, ground_size, 20));
Z_ground = zeros(size(X_ground));
surf(ax_sat_3d, X_ground, Y_ground, Z_ground, 'FaceColor', [0.9, 0.9, 0.9], ...
     'FaceAlpha', 0.3, 'EdgeColor', [0.7, 0.7, 0.7], 'EdgeAlpha', 0.2, ...
     'DisplayName', 'Ground Plane');

% Plot LEO satellite at orbital altitude
plot3(ax_sat_3d, s_x/km_scale, s_y/km_scale, s_z/km_scale, ...
      'bx', 'MarkerSize', 25, 'MarkerFaceColor', 'b', 'DisplayName', ...
      'LEO Satellite');

% Draw vertical line from satellite to ground
plot3(ax_sat_3d, [s_x/km_scale, s_x/km_scale], [s_y/km_scale, s_y/km_scale], ...
      [0, s_z/km_scale], 'b--', 'LineWidth', 1.5, 'HandleVisibility', ...
      'off');

% Add annotation showing altitude
text(ax_sat_3d, s_x/km_scale + 30, s_y/km_scale, s_z/km_scale/2, ...
     sprintf('H = %.0f km', H/km_scale), 'FontSize', 11, 'FontWeight', ...
     'bold', 'Color', 'b');

xlabel(ax_sat_3d, 'X (km)'); ylabel(ax_sat_3d, 'Y (km)'); zlabel(ax_sat_3d, ...
      'Altitude Z (km)');
title(ax_sat_3d, sprintf('LEO Satellite 3D View'));
legend(ax_sat_3d, 'Location', 'northeast');
view(ax_sat_3d, -37.5, 30);

% Use manual axis limits instead of axis equal to show Z clearly
xlim(ax_sat_3d, [-ground_size, ground_size]); 
ylim(ax_sat_3d, [-ground_size, ground_size]); 
zlim(ax_sat_3d, [0, H/km_scale + 50]);

% Adjust data aspect ratio to make Z axis more prominent
% Makes Z axis appear taller relative to X and Y
pbaspect(ax_sat_3d, [1 1 0.4]);

hold(ax_sat_3d, 'off');

% Create text panel for Satellite information
satInfoPanel = uipanel(homeTab, 'Title', 'LEO Satellite Information', ...
                       'Position', [700, 10, 680, 420]);
satInfoText = uitextarea(satInfoPanel, 'Position', [20, 20, 640, 380]);
satInfoText.Editable = 'off';
satInfoText.FontSize = 11;

% Build Satellite information string
satInfoString = sprintf(['Location:\n' ...
    '  ▬ Position: (%.1f, %.1f, %.1f) km\n' ...
    '  ▬ Orbital Altitude: %.0f km\n\n' ...
    ...
    'Technical Specifications:\n' ...
    '  ▬ Channel Bandwidth: %.1f MHz\n' ...
    '  ▬ Transmission Power: %.2f W (%.1f dBm)\n' ...
    '  ▬ Noise Figure: %.1f dB\n' ...
    '  ▬ Path Loss Exponent (η): %.1f\n\n' ...
    ...
    'Network Capacity:\n' ...
    '  ▬ Maximum Network Capacity: %d bbu\n' ...
    '  ▬ Admission Threshold (New Calls): %d bbu\n' ...
    '  ▬ Reserved Capacity (Handoff Calls): %d bbu\n\n' ...
    ...
    'Pricing Parameters:\n' ...
    '  ▬ Base Price (α): %.2f\n' ...
    '  ▬ Network Load Sensitivity (β): %.2f\n' ...
    '  ▬ Available Capacity Sensitivity (γ): %.3f'], ...
    s_x / km_scale, s_y / km_scale, s_z / km_scale, H / km_scale, ...
    W_Sat / 10^6, P_Sat_tx, 10 * log10(P_Sat_tx * 1000), N_Sat, eta, ...
    C_Sat, t_m_Sat, C_Sat - t_m_Sat, ...
    alpha_Sat, beta_Sat, gamma_Sat);


satInfoText.Value = satInfoString;

%% Admitted Users Tab

usersTab = uitab(tabGroup, 'Title', 'Admitted Users');

% Create panel for user list
userListPanel = uipanel(usersTab, 'Title', 'Admitted Users', ...
                        'Position', [10, 10, 250, 820]);

% Find admitted users (not blocked at k=1)
admitted_users = find(all_decisions(:,1) > 0);
num_admitted = length(admitted_users);

% Create listbox for admitted users
userListBox = uilistbox(userListPanel, 'Position', [10, 10, 230, 780]);
userListBox.Items = arrayfun(@(x) sprintf('User %d', x), admitted_users, ...
                            'UniformOutput', false);
userListBox.ItemsData = admitted_users;

% Create panel for plots
userPlotsPanel = uipanel(usersTab, 'Title', 'User Details', ...
                         'Position', [270, 10, 1100, 820]);

% Create axes for network selection plot
ax_user_network = uiaxes(userPlotsPanel, 'Position', [50, 430, 1000, 350]);
title(ax_user_network, 'Select a user to view details');

% Create axes for cost plot
ax_user_cost = uiaxes(userPlotsPanel, 'Position', [50, 30, 1000, 350]);

% Callback for user selection
userListBox.ValueChangedFcn = @(src, event) updateUserPlots(src.Value, ...
              ax_user_network, ax_user_cost, all_decisions, all_costs, K);

%% Network Revenue Tab

revenueTab = uitab(tabGroup, 'Title', 'Network Revenue');

% Create panel for revenue per position update
revPerUpdatePanel = uipanel(revenueTab, 'Position', ...
                           [10, 420, 1360, 410]);
ax_rev_per_update = uiaxes(revPerUpdatePanel, 'Position', ...
                           [50, 30, 1260, 340]);

% Plot Revenue per Update
plot(ax_rev_per_update, 1:K, R_over_time, 'k-x', 'LineWidth', 2, ...
    'MarkerSize', 6, 'MarkerFaceColor', 'k');
xlabel(ax_rev_per_update, 'Position Update Index (k)');
ylabel(ax_rev_per_update, 'Operator Revenue');
title(ax_rev_per_update, 'Network Operator Revenue per Position Update');
grid(ax_rev_per_update, 'on');

% Add labels to data points
for i = 1:K
    text(ax_rev_per_update, i, R_over_time(i), sprintf('%.2f', ...
         R_over_time(i)), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', ...
        'FontSize', 10, 'Color', 'k');
end

% Create panel for cumulative revenue
revCumulativePanel = uipanel(revenueTab, 'Position', ...
                            [10, 10, 1360, 410]);
ax_rev_cumulative = uiaxes(revCumulativePanel, 'Position', ...
                          [50, 30, 1260, 340]);

% Plot Cumulative Revenue
plot(ax_rev_cumulative, 1:K, R_cumulative, 'k-x', 'LineWidth', 2, ...
    'MarkerSize', 6, 'MarkerFaceColor', 'k');
xlabel(ax_rev_cumulative, 'Position Update Index (k)');
ylabel(ax_rev_cumulative, 'Cumulative Revenue');
title(ax_rev_cumulative, 'Cumulative Network Operator Revenue');
grid(ax_rev_cumulative, 'on');

% Add labels to data points
for i = 1:K
    text(ax_rev_cumulative, i, R_cumulative(i), sprintf('%.2f', ...
        R_cumulative(i)), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', ...
        'FontSize', 10, 'Color', 'k');
end

%% QoS Performance Tab

qosTab = uitab(tabGroup, 'Title', 'QoS Performance');

% Create panel for blocking probability
qosBlockPanel = uipanel(qosTab, 'Position', [10, 420, 680, 410]);
ax_qos_block = uiaxes(qosBlockPanel, 'Position', [50, 30, 600, 340]);
stairs(ax_qos_block, 1:K, P_B_over_k, 'r-', 'LineWidth', 2);
xlabel(ax_qos_block, 'Position Update Index (k)');
ylabel(ax_qos_block, 'Blocking Probability P_B');
title(ax_qos_block, 'Call Blocking Probability');
grid(ax_qos_block, 'on');
ylim(ax_qos_block, [0, max(P_B_over_k) + 0.1]);

% Create panel for dropping probability
qosDropPanel = uipanel(qosTab, 'Position', [700, 420, 680, 410]);
ax_qos_drop = uiaxes(qosDropPanel, 'Position', [50, 30, 600, 340]);
stairs(ax_qos_drop, 1:K, P_D_over_k, 'b-', 'LineWidth', 2);
xlabel(ax_qos_drop, 'Position Update Index (k)');
ylabel(ax_qos_drop, 'Dropping Probability P_D');
title(ax_qos_drop, 'Call Dropping Probability');
grid(ax_qos_drop, 'on');
ylim(ax_qos_drop, [0, max(P_D_over_k) + 0.1]);

% % Create panel for Terrestrial Markov Model
% qosTerPanel = uipanel(qosTab, 'Position', [10, 10, 680, 410]);
% ax_qos_ter = uiaxes(qosTerPanel, 'Position', [50, 30, 600, 340]);
% plot(ax_qos_ter, lambda_m_vals, P_B_Ter, 'r-o', 'LineWidth', 2, ...
%     'MarkerSize', 4);
% hold(ax_qos_ter, 'on');
% plot(ax_qos_ter, lambda_m_vals, P_D_Ter, 'b-s', 'LineWidth', 2, ...
%     'MarkerSize', 4);
% hold(ax_qos_ter, 'off');
% xlabel(ax_qos_ter, 'New Call Arrival Rate \lambda_m');
% ylabel(ax_qos_ter, 'Probability');
% title(ax_qos_ter, 'QoS: Terrestrial Network (Markov Model)');
% legend(ax_qos_ter, 'P_B (Blocking)', 'P_D (Dropping)', 'Location', 'best');
% grid(ax_qos_ter, 'on');
% 
% % Create panel for Satellite Markov Model
% qosSatPanel = uipanel(qosTab, 'Position', [700, 10, 680, 410]);
% ax_qos_sat = uiaxes(qosSatPanel, 'Position', [50, 30, 600, 340]);
% plot(ax_qos_sat, lambda_m_vals, P_B_Sat, 'r-o', 'LineWidth', 2, ...
%     'MarkerSize', 4);
% hold(ax_qos_sat, 'on');
% plot(ax_qos_sat, lambda_m_vals, P_D_Sat, 'b-s', 'LineWidth', 2, ...
%     'MarkerSize', 4);
% hold(ax_qos_sat, 'off');
% xlabel(ax_qos_sat, 'New Call Arrival Rate \lambda_m');
% ylabel(ax_qos_sat, 'Probability');
% title(ax_qos_sat, 'QoS: Satellite Network (Markov Model)');
% legend(ax_qos_sat, 'P_B (Blocking)', 'P_D (Dropping)', 'Location', 'best');
% grid(ax_qos_sat, 'on');

%% Helper Functions
%% Check Admission
function [status, L_n, C_n_avail] = check_admission(L_n, C_n_avail, b_k, ...
         C_n, t_m_n, call_type)
    % Projected network load
    L_proj = L_n + b_k;
    
    % New call: Check if the projected load <= new calls admission threshold
    if strcmp(call_type, 'new')
        if L_proj <= t_m_n
            % Mark the network as feasible
            status = 'feasible';
            % Update network state
            L_n = L_proj; C_n_avail = 1 - (L_proj / C_n);
        else
            % Network not feasible → block the call
            status = 'blocked';
        end
    end

    % Handoff call: Check if the projected load <= total network capacity
    if strcmp(call_type, 'handoff')       
        if L_proj <= C_n
            % Mark network as feasible
            status = 'feasible';
            % Update network state
            L_n = L_proj; C_n_avail = 1 - (L_proj / C_n);
        else
            % Network not feasible → drop the call
            status = 'dropped';
        end
    end
end

%% Update User Plots Callback
function updateUserPlots(user_id, ax_network, ax_cost, all_decisions, ...
         all_costs, K)
    % Update network selection plot
    cla(ax_network);
    user_decisions = all_decisions(user_id, :);
    
    % Access coverage status from base workspace
    all_coverage = evalin('base', 'all_coverage');
    user_coverage = all_coverage(user_id, :);
    
    % Create color-coded plot with coverage status background
    colors = zeros(K, 3);
    for k = 1:K
        switch user_decisions(k)
            case 1
                colors(k, :) = [0, 0.5, 0]; % Green for Terrestrial
            case 2
                colors(k, :) = [0, 0, 1]; % Blue for Satellite
            case 0
                colors(k, :) = [1, 0, 0]; % Red for Blocked
            case -1
                colors(k, :) = [0.5, 0, 0.5]; % Purple for Dropped
        end
    end
    
    hold(ax_network, 'on');
    
    % Plot coverage status as background patches
    y_range = [-1.5, 2.5];
    for k = 1:K
        if user_coverage(k) == 1
            % Light green for terrestrial coverage
            patch(ax_network, [k-0.5, k+0.5, k+0.5, k-0.5], ...
                  [y_range(1), y_range(1), y_range(2), y_range(2)], ...
                  [0.8, 1, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        else
            % Light red for satellite-only coverage
            patch(ax_network, [k-0.5, k+0.5, k+0.5, k-0.5], ...
                  [y_range(1), y_range(1), y_range(2), y_range(2)], ...
                  [1, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        end
    end
    
    % Plot network selection using stairs plot
    stairs(ax_network, 1:K, user_decisions, 'k-', 'LineWidth', 2);
    
    % Add markers on the data points
    for k = 1:K
        plot(ax_network, k, user_decisions(k), 'x', 'MarkerSize', 10, ...
             'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', ...
             'LineWidth', 1.5);
    end
    hold(ax_network, 'off');
    
    xlabel(ax_network, 'Position Update Index (k)');
    ylabel(ax_network, 'Network Selection');
    title(ax_network, sprintf('User %d - Network Selection History', ...
          user_id));
    yticks(ax_network, [-1, 0, 1, 2]);
    yticklabels(ax_network, {'Dropped', 'Blocked', 'Terrestrial', ...
               'Satellite'});
    grid(ax_network, 'on');
    ylim(ax_network, y_range);
    
    % Update cost plot - Line graph with total cost only
    cla(ax_cost);
    
    % Get user data from workspace
    user_costs_total = all_costs(user_id, :);
    
    hold(ax_cost, 'on');
    
    % Plot coverage status as background patches
    y_max = max(user_costs_total) * 1.1;
    if y_max == 0
        y_max = 1;
    end
    
    for k = 1:K
        if user_coverage(k) == 1
            % Light green for terrestrial coverage
            patch(ax_cost, [k-0.5, k+0.5, k+0.5, k-0.5], ...
                  [0, 0, y_max, y_max], ...
                  [0.8, 1, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        else
            % Light red for satellite-only coverage
            patch(ax_cost, [k-0.5, k+0.5, k+0.5, k-0.5], ...
                  [0, 0, y_max, y_max], ...
                  [1, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        end
    end
    
    % Plot total cost line
    plot(ax_cost, 1:K, user_costs_total, 'k-x', 'LineWidth', 2, ...
        'MarkerSize', 6);
    
    % Add cost labels on each point (two decimal places)
    y_offset = y_max * 0.03; % 3% of y-axis range for offset
    for k = 1:K
        if user_costs_total(k) > 0
            text(ax_cost, k, user_costs_total(k) + y_offset, ...
                sprintf('%.2f', user_costs_total(k)), 'VerticalAlignment', ...
                'bottom', 'HorizontalAlignment', 'center', ...
                'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
        end
    end
    
    hold(ax_cost, 'off');
    
    xlabel(ax_cost, 'Position Update Index (k)');
    ylabel(ax_cost, 'Total Cost');
    title(ax_cost, sprintf('User %d - Total Cost per Position Update', ...
          user_id));
    grid(ax_cost, 'on');
    ylim(ax_cost, [0, y_max]);
end