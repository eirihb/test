% TTK4135 - Helicopter lab
% Hints/template for problem 2.
% Updated spring 2018, Andreas L. Fl�ten

%% Initialization and model definition
init02; % Change this to the init file corresponding to your helicopter

% Discrete time system model. x = [lambda r p p_dot]'
delta_t	= 0.25; % sampling time
A1 = [[0 1 0            0         0         0];
      [0 0 -K_2         0         0         0];
      [0 0 0            1         0         0];
      [0 0 -K_1*K_pp    -K_1*K_pd 0         0];
      [0 0 0            0         0         1];
      [0 0 0            0         (-K_ep*K_3) (-K_3*K_ed)];] * delta_t + eye(6);
  
B1 = [0        0;
      0        0;
      0        0;
      K_1*K_pp 0;
      0        0;
      0        K_3*K_ep] * delta_t; 

% Number of states and inputs
mx = size(A1,2); % Number of states (number of columns in A)
mu = size(B1,2); % Number of inputs(number of columns in B)

% Initial values
x1_0 = pi;                               % Lambda
x2_0 = 0;                               % r
x3_0 = 0;                               % p
x4_0 = 0;                               % p_dot
x5_0 = 0;                               % e
x6_0 = 0;                               % e_dot
x0 = [x1_0 x2_0 x3_0 x4_0 x5_0 x6_0]';  % Initial values

% Time horizon and initialization
N  = 40;                                % Time horizon for states
M  = N;                                 % Time horizon for inputs
z  = zeros(N*mx+M*mu,1);                % Initialize z for the whole horizon
z0 = z;                                 % Initial value for optimization
z0(1) = pi;

% Bounds
ul 	    = -30*pi/180;                   % Lower bound on control pitch
uu 	    =  30*pi/180;                   % Upper bound on control ptich
ul    = ul*ones(mu,1);
uu   = uu*ones(mu,1);


xl      = -Inf*ones(mx,1);              % Lower bound on states (no bound)
xu      = Inf*ones(mx,1);               % Upper bound on states (no bound)
xl(3)   = ul(1);                           % Lower bound on state x3
xu(3)   = uu(1);                           % Upper bound on state x3
xl(5)   = -30*pi/180;                           % Lower bound on state x5
xu(5)   = 40*pi/180;                           % Upper bound on state x5

%zl = repmat([xl ; ulow],N,1);
%zu = repmat([xu ; uhigh],N,1);


% Generate constraints on measurements and inputs

[vlb,vub]       = gen_constraints(N,M,xl,xu,ul,uu); % hint: gen_constraints
vlb(N*mx+M*mu)  = 0;                    % We want the last input to be zero
vub(N*mx+M*mu)  = 0;                    % We want the last input to be zero

% Generate nonlinear constraint
nonlinearConstraint = @nonLinCon;


% Generate objective function
%        lambda, r, p, p_dot, e, e_dot
Q1 = diag([1,     0, 0,  0,    0, 0]);                           % Weight on state x
P1 = [1 0 ; 0 1];                                % Weight on input
Q = gen_q(Q1,P1,N,M);                                  % Generate Q, hint: gen_q


%objFun = @objectiveFunc;
objFun = @(X) X' * Q * X;
%% Generate system matrixes for linear model
Aeq = gen_aeq(A1,B1,N,mx,mu); % Discretized A, Discretized B, Time horizon, #states, #inputs            % Generate A, hint: gen_aeq
beq = zeros(size(B1,1) *N,1);             % Generate b
beq(1:mx) = x0;

%% Solve quadratic problem with nonlinear consterainsts model
maxIterations = 1000*(N^3);
OPTIONS = optimoptions('fmincon', 'Algorithm','sqp','MaxSQPIter',maxIterations);
tic
z = fmincon(objFun,z0,[],[],Aeq,beq,vlb,vub,nonlinearConstraint,OPTIONS); % hint: fmincon. Type 'doc fmincon' for more info 
t1=toc;
disp('Start plotting');

% Calculate objective value
% phi1 = 0.0;
% PhiOut = zeros(N*mx+M*mu,1);
% for i=1:N*mx+M*mu
%   phi1=phi1+Q(i,i)*z(i)*z(i);
%   PhiOut(i) = phi1;
% end


%%% u blir n� innholdet for begge U-ene. Splitt det opp!
%% Extract control inputs and states
u1  = [z(N*mx+1:2:N*mx+M*2);z(N*mx+M*mu-1)]; % Control input from solution
u2  = [z(N*mx+2:2:N*mx+M*2);z(N*mx+M*mu)]; % Control input from solution

x1 = [x0(1);z(1:mx:N*mx)];              % State x1 from solution
x2 = [x0(2);z(2:mx:N*mx)];              % State x2 from solution
x3 = [x0(3);z(3:mx:N*mx)];              % State x3 from solution
x4 = [x0(4);z(4:mx:N*mx)];              % State x4 from solution
x5 = [x0(5);z(5:mx:N*mx)];              % State x5 from solution
x6 = [x0(6);z(6:mx:N*mx)];              % State x6 from solution



num_variables = 8/delta_t;
zero_padding = zeros(num_variables,1);
unit_padding  = ones(num_variables,1);

u1   = [zero_padding; u1; zero_padding];
u2   = [zero_padding; u2; zero_padding];
x1  = [pi*unit_padding; x1; zero_padding];
x2  = [zero_padding; x2; zero_padding];
x3  = [zero_padding; x3; zero_padding];
x4  = [zero_padding; x4; zero_padding];
x5  = [zero_padding; x5; zero_padding];
x6  = [zero_padding; x6; zero_padding];


%% Plotting
t = 0:delta_t:delta_t*(length(u1)-1);

figure(2)
subplot(811)
stairs(t,u1),grid
ylabel('u1')
subplot(812)
stairs(t,u2),grid
ylabel('u2')
subplot(813)
plot(t,x1,'m',t,x1,'mo'),grid
ylabel('lambda')
subplot(814)
plot(t,x2,'m',t,x2','mo'),grid
ylabel('r')
subplot(815)
plot(t,x3,'m',t,x3,'mo'),grid
ylabel('p')
subplot(816)
plot(t,x4,'m',t,x4','mo'),grid
xlabel('tid (s)'),ylabel('pdot')
subplot(817)
plot(t,x5,'m',t,x5','mo'),grid
xlabel('tid (s)'),ylabel('e')
subplot(818)
plot(t,x6,'m',t,x6','mo'),grid
xlabel('tid (s)'),ylabel('edot')

trajectoryu = [t' u1 u2];

trajectoryx = [t' x1 x2 x3 x4 x5 x6];


%% LQR
%           lambda, r, p, pdot, e, edot
Q_lqr = diag([10,.1,100,0.1, 10, 0.1]);
%             p_c, e_c
R_lqr = diag([10 , 1]);

K_lqr = dlqr(A1, B1, Q_lqr, R_lqr);

% 
% %% LQR
% %           lambda, r, p, pdot, e, edot
% Q_lqr = diag([0.000000000001,0.0000000000000000001, 0.01,1000, 0.00000000000000001, 0.00000000000000000000000000000000000001]);
% %             p_c, e_c
% R_lqr = diag([0.000000000000000001 , 1]);
% 
% K_lqr = dlqr(A1, B1, Q_lqr, R_lqr);







