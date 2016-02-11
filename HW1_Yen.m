%% Numerical Method HW1
%% Bing-Jie Yen
clear all;
clc;
%% Discrete-time version of a Mortensen-Pissarides with aggregate fluctuation
%% 1. Set parameters (global variables)
global delta alpha A rho_z sig_epsilon mu b kappa beta z N sig_z
delta =0.081;           % Seperation rate
alpha =0.72;            % Elasticity of matching
A     = 0.158;          % Mactching efficientcy
rho_z =0.9895;          % Autocorrelation of weekly productivity
sig_epsilon = 0.0034; % Standard deviation for innovations
mu    = 0.72;           % Bargaining weight for workers
b     = 0.4;            % Unempoyment utility
kappa = 0.6;            % Posting cost is about 3 days of output
beta  = 0.999;          % Weekly discount rate



%% 2. Functional forms (Matching function)
%m_fun = min(1, A*(theta.^(1-alpha)));

%% Where the theta comes from? 
%% 3. Solving Mortensen-Pissarides model 
%% 4. Tauchen's method (follow the lecture notes step by step...)
% Goal: Choosing values for the realization :z and the transition matrix so the markov chain closely mimics AR(1)

% (1) Let rho and sig_epsilon represent autoregressive parameter and
% standard deviation of innovations : rho sig_epsilon
% (2)Choose the number of grids and set up the relationship between
% variances : N
% (3)Select the size of the grid by choosing lambda :dz(i)

lambda = 3;
sig_z= sig_epsilon/(sqrt(1-(rho_z)^2));
z_min= -lambda.*sig_z;  
z_max = lambda.*sig_z;
N=50;        
% (4) Creating equidistant z

step= (-2.*z_min)/(N-1);
z=zeros(N,1);
for i=1:N 
    z(i)= z_min +(i-1).*step;
end


% (5) In Tauchen's method, let dz denote the half of the distance between two consecutive grid points


%dz= zeros(N-1,1);  % Snce N points include N-1 midpoints
%for i= 1:(N-1)
%    dz(i)= (z(i+1)+z(i))/2;
%end
%% 

%% (5) Compute the transition matrix p_{ij}  
P=zeros(N,N);
for i= 1:N
    P(i,1)= normcdf((z_min-rho_z.*z(i))/sig_epsilon+step/(2.*sig_epsilon));
    
        for j=2:N-1
            P(i,j)= normcdf((z(j)-rho_z.*z(i))/sig_epsilon+step/(2.*sig_epsilon)- normcdf((z(j)-rho_z.*z(i))/sig_epsilon-step/(2.*sig_epsilon));
        end
        for j=1:N-1
            P(i,N)=1-sum(P(i,j));
        end          
end
global P
%% 5. Solve the market tightness from MP model ( which is in MP_model.m )
%% 
MP_root=@MP_model;
theta_guess=zeros(N,1);
theta =fsolve(MP_root,theta_guess);
global theta 

% Find market tightness: theta (1) Matching Function (2)MP-function
global p_theta q_theta
p_theta=MatchingFunction(theta); 
q_theta=p_theta./theta;



end 

















% 3. Tauchen's method
% 4. 