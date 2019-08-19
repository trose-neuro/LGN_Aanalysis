function [pval chi2stat tbl] = propstat(n1, N1, n2, N2)
% TR 2015 -> 
% https://de.mathworks.com/matlabcentral/answers/96572-how-can-i-perform-a-chi-square-test-to-determine-how-statistically-different-two-proportions-are-in
% http://www.evernote.com/l/AATcdR8_-AxPf4s5isQPAtcV1IzP_Cdr1Yc/

x1 = [repmat('a',N1,1); repmat('b',N2,1)];

x2 = [repmat(1,n1,1);
    
repmat(2,N1-n1,1);

repmat(1,n2,1);

repmat(2,N2-n2,1)];

[tbl,chi2stat,pval] = crosstab(x1,x2);
