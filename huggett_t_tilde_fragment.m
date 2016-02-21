clear all
close all
% PROGRAM NAME: ps4huggett.m
clear, clc

% PARAMETERS
beta = .9932; %discount factor 
sigma = 1.5; % coefficient of risk aversion
b = 0.5; % replacement ratio (unemployment benefits)
y_s = [1, b]; % endowment in employment states
PI = [.97 .03; .5 .5]; % transition matrix


% ASSET VECTOR
a_lo = -2; %lower bound of grid points
a_hi = 3; %upper bound of grid points
num_a = 700;

a = linspace(a_lo, a_hi, num_a); % asset (row) vector

% INITIAL GUESS FOR q
q_min = 0.98;
q_max = 1.1;
%num_q = 10;
%q = linspace(q_min, q_max, num_q);
q_guess = (q_min + q_max) / 2;

% ITERATE OVER ASSET PRICES
aggsav = 1 ;
i=1;
while abs(aggsav) >= 0.01 ;
    
    q_guess = (q_min + q_max) / 2;
    
    % CURRENT RETURN (UTILITY) FUNCTION
    cons = bsxfun(@minus, a', q_guess * a);
    cons = bsxfun(@plus, cons, permute(y_s, [1 3 2]));
    ret = (cons .^ (1-sigma)) ./ (1 - sigma); % current period utility
    ret(cons < 0) = -Inf;
    
    % INITIAL VALUE FUNCTION GUESS
    v_guess = zeros(2, num_a);
    
    % VALUE FUNCTION ITERATION
    v_tol = 1;
    while v_tol >.0000001;
        % CONSTRUCT TOTAL RETURN FUNCTION
        value_mat = ret + beta * ...
        repmat(permute((PI * v_guess), [3 2 1]), [num_a 1 1]);
        % CHOOSE HIGHEST VALUE (ASSOCIATED WITH a' CHOICE)
        [vfn, pol_indx] = max(value_mat, [], 2);
        vfn =permute(vfn, [3,1,2]);
        
        % what is the distance between current guess and value function
        v_tol = max(max(abs(vfn - v_guess)));
    
        % if distance is larger than tolerance, update current guess and
        % continue, otherwise exit the loop
        v_guess = vfn;
  
    end;
    
    % KEEP DECSISION RULE
    pol_indx = permute(pol_indx, [3, 1, 2]);
    pol_fn = a(pol_indx);
   
    % SET UP INITITAL DISTRIBUTION
    Mu=ones(2, num_a)/(2*num_a);
    
    mu_tol = 1;
    while mu_tol >.0000001;
    % ITERATE OVER DISTRIBUTIONS
    [emp_ind, a_ind, mass] =find(Mu); % find non-zero indices
    
    MuNew = zeros(size(Mu));
    for ii = 1:length(emp_ind)
        apr_ind = pol_indx(emp_ind(ii), a_ind(ii));
        MuNew(:, apr_ind) = MuNew(:, apr_ind) + ...
            (PI(emp_ind(ii), :) * mass(ii))';
        
    end
    
    mu_tol = max(max(abs(MuNew - Mu)));
    
    Mu=MuNew;
    
    end
    
    aggsav = sum(sum(Mu.*pol_fn))
    
    
    if aggsav>0
        q_min=q_guess;
    else
        q_max=q_guess;
    end
    
    
    i=i+1
        
end


% Plot the histogram as the note

wealth=repmat(a, [2,1])+repmat(y_s', [1, num_a]);
plot_e=bar(wealth(1,:), Mu(1,:), 'FaceColor', 'b', 'EdgeColor', 'b', 'LineWidth', 0.1)
hold on
plot_u=bar(wealth(2,:), Mu(2,:), 'FaceColor', 'r', 'EdgeColor', 'r', 'LineWidth', 0.1)
hold off
title('Wealth distribution of a=-2 and \pi(u|u)=0.5')
xlabel('Wealth')
ylabel('Fraction of the population')
Legend('Employed', 'Unemployed')

% Plot Lorenz curves
figure
wealth_total=[wealth(1,:), wealth(2,:)]
Mu_total=[Mu(1,:), Mu(2,:)]
[wealth_sort wealth_indx]=sort(wealth_total)
Mu_sort=Mu(wealth_indx)
wealth_total_sort=wealth_sort.*Mu_sort
Mu_sort_sum=cumsum(Mu_sort)
wealth_sum=sum(wealth_total_sort)
wealth_total_sort=cumsum(wealth_total_sort)
wealth_total_sort=wealth_total_sort./wealth_sum
plot(Mu_sort_sum, wealth_total_sort)
title('Lorenz curve of wealth')
xlabel('Cumulative shares of people from lowest to highest wealth')
ylabel('Cumulative shares of wealth')


figure
earning=repmat(y_s', [1, num_a])
earning=[earning(2,:), earning(1,:)]
Mu_earning=[Mu(2,:), Mu(1,:)]
earning_total_sort=earning.*Mu_earning
Mu_earning_sum=cumsum(Mu_earning)
earning_sum=sum(earning_total_sort)
earning_total_sort=cumsum(earning_total_sort)
earning_total_sort=earning_total_sort./earning_sum
plot(Mu_earning_sum, earning_total_sort)
title('Lorenz curve of earnings')
xlabel('Cumulative shares of people from lowest to highest earnings')
ylabel('Cumulative shares of earnings')


% Gini coefficient of wealth is 0.62. Gini coefficient of earning is 0.03
ope=[0:1:2*num_a-1]

wealth_total_sort=wealth_sort.*Mu_sort
wealth_total_sort=cumsum(wealth_total_sort)
wealth_total_sort_1=zeros(1,2*num_a)
wealth_total_sort_1(1,2:2*num_a)=wealth_total_sort(ope(1,2:2*num_a))
wealth_total_sort_g=wealth_total_sort_1+wealth_total_sort
wealth_total_sort_g=sum(wealth_total_sort_g.*Mu_sort)
G_wealth=1-(wealth_total_sort_g/wealth_sum)


earning_total_sort=earning.*Mu_earning
earning_total_sort=cumsum(earning_total_sort)
earning_total_sort_1=zeros(1,2*num_a)
earning_total_sort_1(1,2:2*num_a)=earning_total_sort(ope(1,2:2*num_a))
earning_total_sort_g=earning_total_sort_1+earning_total_sort
earning_total_sort_g=sum(earning_total_sort_g.*Mu_earning)
G_earning=1-(earning_total_sort_g/earning_sum)


% Question 4
% Since long run pi_u is 0.0566, the consumption with insurance is 0.9717
c_ins=0.9434+.5*.0566;

% The life time utility of the consumption with insurance W_fb is -298.37
ret_ins=(1/(1-beta))*(((c_ins)^(1-sigma))/(1-sigma));

% The calculation of lambda
lambda=((vfn.^(-1)).*ret_ins).^(1/(1-sigma))-1;

% The fraction of household benefit from the plan is 54%
welf_inc=(lambda>0);
welf_inc_per=sum(sum(welf_inc.*Mu));

% The economy-wide welfare gain 0.001353
welf=sum(sum(lambda.*Mu));