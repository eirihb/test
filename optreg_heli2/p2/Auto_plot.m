filename1 = 'p2p4_1.mat';
m1 = matfile(filename1);
array = m1.ulambdarppdot;

%%Time
time = array(1,:);

%Measured vals
pitch_setpoint = array(2,:);
lambda = array(3,:);
lambda_dot = array(4,:);
pitch = array(5,:);
pitch_rate = array(6,:);
figure;

%%u
subplot(5,1,1);
hold on;
xlim([0 15]);
ylim([-6 11]);
plot(time, pitch_setpoint, 'r');
%handle=legend('u');
%set(handle,'Interpreter','latex')
xlabel('Time t [s]');
ylabel('u [V]');

%%lambda
subplot(5,1,2);
hold on;
xlim([0 15]);
ylim([-2 1]);
plot(time, lambda, 'r');
%legend('{\lambda}');
xlabel('Time t [s]');
ylabel('{\lambda} [1]');

%%Lambda_dot
subplot(5,1,3);
hold on;
xlim([0 15]);
plot(time, lambda_dot, 'r');
%hl = legend('$\dot{\lambda}$');
%set(hl,'Interpreter','latex')
xlabel('Time t [s]');
ylabel('$\dot{\lambda} [s^{-1}]$','Interpreter','latex');

%%Pitch
subplot(5,1,4);
hold on;
xlim([0 15]);
plot(time, pitch, 'r');
%hl1=legend('$p$');
%set(hl1,'Interpreter','latex')
xlabel('Time t [s]');
ylabel('$p [1]$','Interpreter','latex');

%%Pitch rate
subplot(5,1,5);
hold on;
xlim([0 15]);
plot(time, pitch_rate, 'r');
%hl2=legend('$\dot{p}$');
%set(hl2,'Interpreter','latex')
xlabel('Time t [s]');
ylabel('$\dot{p} [s^{-1}]$','Interpreter','latex');
