function [s_n, i_n, r_n, v_n, d_n] = action_sirvd(s, i, r, v, d, beta, gamma, alpha, sigma, p)
% Advance an SIR model one timestep
%
% Usage
%   [s_n, i_n, r_n] = action_sir(s, i, r, beta, gamma)
% 
% Arguments
%   s = current number of susceptible individuals
%   i = current number of infected individuals
%   r = current number of recovered individuals
%   v = current number of vaccinated individuals
%   d = current number of dead individuals
%   
%   beta = infection rate parameter
%   gamma = recovery rate paramter
%   alpha = reinfection rate
%   sigma = vaccination rate
%   p = percent dying from infection
% 
% Returns
%   s_n = next number of susceptible individuals
%   i_n = next number of infected individuals
%   r_n = next number of recovered individuals
%   v_n = next number of vaccinated individuals
%   d_n = next number of deceased individuals

% compute new infections and recoveries
infected = beta * i * s;
resuseptible = alpha * r;
vaccinated = sigma * s;
recovered = (1-p) * gamma * i;
dead = p * gamma * i;

infected = min(infected, s);
recovered = min(recovered, i);
resuseptible = min(resuseptible, r);

% Update state
s_n = s - infected + resuseptible - vaccinated;
i_n = i + infected - recovered - dead;
r_n = r + recovered - resuseptible;
v_n = v + vaccinated;
d_n = d + dead;

% Enforce invariants; necessary since we're doing a discrete approx.
s_n = max(s_n, 0);
i_n = max(i_n, 0);
r_n = max(r_n, 0);
    
end