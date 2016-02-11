function p_theta=MatchingFunction(theta)
global alpha A theta 

p_theta=min(1,A*(theta.^(1-alpha)));
end