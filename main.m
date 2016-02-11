%% Numerical Method HW1
%% Bing-Jie Yen
clear all;
clc;
%% Discrete-time version of a Mortensen-Pissarides with aggregate fluctuation
%% 1. Set parameters (global variables)%% Discrete-time version of a Mortensen-Pissarides with aggregate fluctuation
%% 1. Set parameters (global variables)
%% 2. Space discrete AR(1) : Tauchen method
%% 3. Find the roots by solving stochastic difference equation  

global delta alpha A rho_z sig_epsilon mu b kappa beta z N lambda
delta =0.081;           % Seperation rate
alpha =0.72;            % Elasticity of matching
A     = 0.158;          % Mactching efficientcy
rho_z =0.9895;          % Autocorrelation of weekly productivity
sig_epsilon = 0.0034; % Standard deviation for innovations
mu    = 0.72;           % Bargaining weight for workers
b     = 0.4;            % Unempoyment utility
kappa = 0.6;            % Posting cost is about 3 days of output
beta  = 0.999;          % Weekly discount rate
N =50;                   % Number of grids
lambda = 3;


%% 2. Tauchen's method (follow the lecture notes step by step...)
% Goal: Choosing values for the realization :z and the transition matrix so the markov chain closely mimics AR(1)
% (1) Let rho_z and sig_epsilon represent autoregressive parameter and standard deviation of innovations : rho sig_epsilon
% (2)Choose the number of grids and set up the relationship between variances : N
% (3)Select the size of the grid by choosing lambda :dz(i)
% (4) Creating equidistant z
% (5) Compute the transition matrix p_{ij}  

sig_z = sig_epsilon/(sqrt(1-(rho_z)^2));
z_min= -lambda.*sig_z;  
z_max = lambda.*sig_z;


step= (-2.*z_min)/(N-1);                                             % (4) Creating equidistant z
z=zeros(N,1);
for i=1:N 
    z(i)= z_min +(i-1).*step;
end


% (5) Compute the transition matrix p_{ij}  
P =zeros(N,N);
for i= 1:N
    P(i,1)= normcdf((z_min-rho_z.*z(i))/sig_epsilon+step/(2.*sig_epsilon));

        for j=2:N-1
            P(i,j)= normcdf((z(j)-rho_z.*z(i))/sig_epsilon + step/(2.*sig_epsilon))- normcdf((z(j)-rho_z.*z(i))/sig_epsilon-step/(2.*sig_epsilon));
        end

        for j=1:N-1
            P(i,N)=1-sum(P(i,j));
        end          
end


%% 3. Solve the market tightness 
%(1) Refer: http://www.mathworks.com/help/optim/ug/fsolve.html


MP_root= @(theta) P*(1-mu)*(z-b)-kappa*mu*exp(x)+(kappa/A)*(1-delta)*exp(alpha*x)-(kappa/A*beta)*exp(alpha*x)
log_theta=zeros(N,1);
problem.options = optimoptions('fsolve','Display','none','PlotFcns',@optimplotfirstorderopt);
problem.objective = MP_root;
problem.x0 = log_theta;
problem.solver = 'fsolve';

log_theta_replace=fsolve(problem);


figure;
subplot(2,1,1);
plot(log_theta_replace);
xtitle('Solution to Market Tightness- log(theta)');
xlabel('the # of grid');
ylabel('log(theta)');

subplot(2,1,2)
plot(MP_root(log_theta_replace));
title('Residual');
xlabel('the # of grid');
ylabel('residual');

%

%% 4. Simulate from the Markov Chain

% (2) The revised question : solving log(z) instead of z
%  log Z_{t+1} = rho log Z_t + epsilon
%First, define a function for log Z

log_mean(Z)= 0; % the productivity mean is 1, therefore, log (mean z)=0
logZ= 


%%
end 























% 3. Tauchen's method
% 4. 
