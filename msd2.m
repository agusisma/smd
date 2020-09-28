%% Mass-Spring-Damper System
% By Agus Hasan

clear;
clc;

dt = 0.001;                 % time step
tf = 20;                    % simulation time

m = 10;                     % mass
k = 500;                    % spring coefficient
b = 5;                      % damper coefficient

% system description
B = [0;1/m]*dt;             % control matrix
H = [1 0 0 0;               % measurement matrix digital twin
     0 1 0 0];
C = [1 0;                   % measurement matrix physical twin
     0 1];
u = 0;                      % initial input

% error covariance matrix
Q = 0.1*eye(4);             % model covariance matrix
R = 0.1;                    % measurement covariance matrix

% initial data
x     = [10 0]';            % initial condition
xhat  = [10 0 0 400]';      % estimated initial condition
Pplus = 10000*eye(4);       % initial matrix propagation

% for plotting
xArray    = [];
xhatArray = [];

for i=1:tf/dt
    
    if i>5000
        b = 6;
    end
    if i>10000
        u = 50000;
    end
    if i>10010
        u = 0;
    end    
    xArray    = [xArray x];
    xhatArray = [xhatArray xhat];
    % Simulate the system
    x = [1 dt;-k*dt/m 1-(b*dt/m)]*x+B*u;
    y = C*x;
    % Prediction
    F = [1 dt 0 0;          % Jacobian matrix
         -xhat(4)*dt/m 1-(xhat(3)*dt/m) -dt*xhat(2)/m -xhat(1)*dt/m;
         0 0 1 0;
         0 0 0 1];
    xhat  = [xhat(1)+xhat(2)*dt;(-xhat(4)*dt/m)*xhat(1)+(1-(xhat(3)*dt/m))*xhat(2)+dt*u/m;xhat(3);xhat(4)];
    Pmin  = F*Pplus*F' + Q;
    % Update
    K     = Pmin*H'*inv(H*Pmin*H' + R);
    Pplus = (eye(4)-K*H)*Pmin;
    xhat  = xhat + K*(y-H*xhat);
end

figure(1)
subplot(2,1,1)
plot(dt:dt:tf,xArray(1,:),'-b','LineWidth',3)
hold on;
plot(dt:dt:tf,xhatArray(1,:),':r','LineWidth',3)
legend('physical twin position','digital twin position')
set(gca,'FontSize',24)
grid on;
grid minor;
subplot(2,1,2)
plot(dt:dt:tf,xArray(2,:),'-b','LineWidth',3)
hold on;
plot(dt:dt:tf,xhatArray(2,:),':r','LineWidth',3)
legend('physical twin velocity','digital twin velocity')
set(gca,'FontSize',24)
grid on;
grid minor;

figure(2)
subplot(2,1,1)
plot(dt:dt:tf,[5*ones(1,5000) 6*ones(1,15000)],'-b','LineWidth',3)
hold on;
plot(dt:dt:tf,xhatArray(3,:),':r','LineWidth',3)
legend('physical twin damping coef.','digital twin damping coef.')
ylim([4.5 6.5])
set(gca,'FontSize',24)
grid on;
grid minor;
subplot(2,1,2)
plot(dt:dt:tf,[500*ones(1,5000) 500*ones(1,15000)],'-b','LineWidth',3)
hold on;
plot(dt:dt:tf,xhatArray(4,:),':r','LineWidth',3)
ylim([405 550])
legend('physical twin spring coef.','digital twin spring coef.')
set(gca,'FontSize',24)
grid on;
grid minor;

figure(3)
curve1 = animatedline('Color','b','LineWidth',2);
curve2 = animatedline('Color','r','LineStyle',':','LineWidth',2);
set(gca,'XLim',[0 20],'YLim',[-15 15]);
legend('physical twin position','digital twin position')
ylabel('position (m)')
xlabel('time (s)')
grid on;
tm = dt:10*dt:tf;
for i = 1:length(tm)
    addpoints(curve1,tm(i),xArray(1,10*i));
    hold on
    addpoints(curve2,tm(i),xhatArray(1,10*i));
    drawnow
    G(i) = getframe(gcf);
end
video = VideoWriter('MSDvel.');
open(video)
writeVideo(video,G)
close(video)

figure(4)
curve2 = animatedline('Color','r','LineStyle',':','LineWidth',2);
set(gca,'XLim',[0 20],'YLim',[4.5 6.5]);
legend('digital twin damping coef.')
ylabel('damping coeff')
xlabel('time (s)')
grid on;
tm = dt:10*dt:tf;
for i = 1:length(tm)
    addpoints(curve2,tm(i),xhatArray(3,10*i));
    drawnow
    G(i) = getframe(gcf);
end
video = VideoWriter('MSDspring.');
open(video)
writeVideo(video,G)
close(video)
