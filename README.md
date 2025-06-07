# unequal-Alloction
Phase II trial design using unequal allocation ratios



> p0 <- 0.1
> p1 <- 0.3
> alpha_star <- 0.1
> power_star <- 0.9
> gamma <- binAlloc(p0,p1);gamma
[1] 1.527525
> 
> 
> singleStageII(p0,p1,alpha_star,power_star,1)
   p0  p1 pow_star alpha_star gamma m1 n1  ss     Power      Alpha    Theta
1 0.1 0.3      0.9        0.1     1 58 58 116 0.9008774 0.05157159 3.857143

> singleStageII(p0,p1,alpha_star,power_star,gamma)
   p0  p1 pow_star alpha_star    gamma m1 n1  ss     Power      Alpha    Theta
1 0.1 0.3      0.9        0.1 1.527525 42 64 106 0.9043931 0.06109451 3.857143

> randPhaseII(p0,p1,alpha_star,power_star,1)
   p0  p1 pow_star alpha_star gamma n1 m1 n2 m2  ss  Power  Alpha   PET DES
1 0.1 0.3      0.9        0.1     1 27 27 32 32 118 0.9013 0.0569 0.408  MM
2 0.1 0.3      0.9        0.1     1 20 20 41 41 122 0.9002 0.0587 0.393 OPT

> randPhaseII(p0,p1,alpha_star,power_star,gamma)
   p0  p1 pow_star alpha_star    gamma n1 m1 n2 m2  ss  Power  Alpha   PET DES
1 0.1 0.3      0.9        0.1 1.527525 32 21 32 21 106 0.9016 0.0603 0.481  MM
2 0.1 0.3      0.9        0.1 1.527525 26 17 40 26 109 0.9012 0.0593 0.483 OPT




