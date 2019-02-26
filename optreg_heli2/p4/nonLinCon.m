function [c, ceq]=nonLinCon(z)
    mx = 6;
    mu = 2;
    sizeOfTimestep = mu + mx;
    N = size(z,1)/sizeOfTimestep;
    alpha = 0.2;
    beta = 20;
    lambda_t=2*pi/3;
    c=zeros(N);
    for n = 1:N
        %                           lambda(n)                       e(n)
        c(n) = alpha * exp(-beta*(z(1+(n-1)*mx)-lambda_t)^2) - z(5+mx*(n-1));
      

    end    
    ceq =[];
end