function phi=objectiveFunc(z)
disp('Objfunc - start');
    sizeOfTimestep = 8;
    mx=6;
    mu=2;
    q1 = 1;
    q2 = 1;
    N = size(z,1)/sizeOfTimestep;    
    %lambda_f = z(1+(N-1)*mx);
    lambda_f = 0;
    phi = 0;
    for i=1:N
        %            lambda_i
        phi = phi + (z(1+(i-1)*mx)-lambda_f)^2 + ...
            ...%   pc_i = u(1)_i
            q1*(z(1+(i-1)*mu + N*mx))^2 + ...
            ...%   pe_i = u(2)_i
            q2*(z(2+(i-1)*mu + N*mx))^2;

    end
   disp('Objfunc - end');
   phi
end

