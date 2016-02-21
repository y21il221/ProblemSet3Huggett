1. The recursive problem is stated in the code.

2. The equilibrium q is 0.994. Therefore, the interest rate is 0.6% per period. Since in the Huggett's model, there are six periods per year, the annual interest rate is 3.6%. The distribution of wealth of the employed group and the unemployed group as the note are also ploted.

3. The Gini coefficient of wealth is 0.62, and the Gini coefficient of earning is 0.03

4. The code of question 4 is in the m-file,

"% Question 4
% Since long run pi_u is 0.0566, the consumption with insurance is 0.9717
c_ins=0.9434+.5*.0566;

% The life time utility of the consumption with insurance W_fb is -298.37
ret_ins=(1/(1-beta))*(((c_ins)^(1-sigma))/(1-sigma));

% The calculation of lambda
lambda=((vfn.^(-1)).*ret_ins).^(1/(1-sigma))-1;

% The fraction of household benefit from the plan is 54%
welf_inc=(lambda>0);
welf_inc_per=sum(sum(welf_inc.*Mu));

% The economy-wide welfare gain is 0.001353
welf=sum(sum(lambda.*Mu));"
