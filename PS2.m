% Question (1)
% state variable: current capitial and productivity 
% control variable: current consumption and next period's capital

close all
%%%% Set up parameters
alpha = 0.35;
beta = 0.99;
delta = 0.025;
sigma = 2;
a=[1.1, 0.678];
prob=[0.977,1-0.977;1-0.926,0.926]


%%%% Set up discretized state space
k_min = 0;
k_max = 45;
num_k = 1000; % number of points in the grid for k

k = linspace(k_min, k_max, num_k); 
%Create a vector of num_k evenly spaced points in the interval[k_min,k_max].

k_mat = repmat(k', [1 num_k]); % this will be useful in a bit


%%%% Set up consumption and return function
% 1st dim(rows): k today, 2nd dim (cols): k' chosen for tomorrow
cons1 =a(1)*(k_mat .^ alpha) + (1 - delta) * k_mat - k_mat'; 
cons2 =a(2)*(k_mat .^ alpha) + (1 - delta) * k_mat - k_mat'; 

ret1 = cons1 .^ (1 - sigma) / (1 - sigma); % return function
ret2 = cons2 .^ (1 - sigma) / (1 - sigma);

% negative consumption is not possible -> make it irrelevant by assigning
% it very large negative utility
ret1(cons1 < 0) = -Inf;
ret2(cons2 < 0) = -Inf;

%%%% Iteration
dis = 1; tol = 1e-06; % tolerance for stopping 
v_guess = zeros(2, num_k);
while dis > tol
    % compute the utility value for all possible combinations of k and k':
    value_mat1 = ret1 + beta * repmat(prob(1,:)*v_guess, [num_k 1]);
    value_mat2 = ret2 + beta * repmat(prob(2,:)*v_guess, [num_k 1]);
    
    % find the optimal k' for every k:
    [vfn1, pol_indx1] = max(value_mat1, [], 2);
    [vfn2, pol_indx2] = max(value_mat2, [], 2);
    
    vfn=[vfn1';vfn2'];
   
    
    % what is the distance between current guess and value function
    dis = max(abs(vfn - v_guess));
    
     % if distance is larger than tolerance, update current guess and
    % continue, otherwise exit the loop
    v_guess = vfn;
    
    
end

g1=k(pol_indx1); % policy function if A=1.1
g2=k(pol_indx2); % policy function if A=0.678



% Question 2
plot(k,vfn1)
xlabel('k') 
ylabel('V(k)')
title('value function over K for each A')
hold on;

plot(k,vfn2)
legend('A_t=1.1','A_t=0.678');

%%%% Ans: Yes, they are increasing and concave

figure

%Question 3
plot(k,g1)
xlabel('k_t') 
ylabel('k_t+1')
title('policy function over K for each A')
hold on;

plot(k,g2)
legend('A_t=1.1','A_t=0.678');
%%%% Ans: Yes, it is increasing in k and A.

figure

s1=g1-(1-delta)*k;
plot(k,s1)
xlabel('k')
ylabel('saving')
title('saving over K for each A')
hold on; 

s2=g2-(1-delta)*k;
plot(k,s2)
legend('A_t=1.1','A_t=0.678');
%%%% Ans: it is incresing in A, but decreasing in K

% Question 4
% Simulation
T_sim=5000;

% draw random unmber for A
rng(1)

rand_nums=rand(T_sim,1);

% turn random numbers into values for A using transition matrix
A_sim=zeros(T_sim,1);
A_sim(1)=1;

for t=1:T_sim;
  if A_sim(t)==1;
    A_sim(t+1)=1 % if rand_nums(t)<pi_hh, else A_sim(t+1)=2
   else; % A_sim(t)
 A_sim(t+1)=2  
    
    

% start with arbitrary capital stock, then follow the policy function
% according to simulated state in the current period
k_sim_index=zeros(T_sim,1)
k_sim_index(1)=5; % capital for the first period
for t=1;
  k_sim_index(t+1)=pol_indx(A_sim(t), k_sim_index(t));
end
% Then, calculate output
