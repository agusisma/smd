clear;
clc;

dt = 0.001;
tf = 10;

m = 10;
k = 500;
b = 5;

% system description
B = [0;1/m]*dt;
H = [1 0 0;
     0 1 0];
C = [1 0;
     0 1];
u = 0;

% error covariance matrix
Q = 0.1*eye(3);
R = 0.1;

% initial data
x     = [10 0]';
xhat  = [10 0 0]';
Pplus = 10000*eye(3);

% for plotting
xArray    = [];
xhatArray = [];

for i=1:tf/dt
    
    if i>5000
        b = 6;
    end
    
    xArray    = [xArray x];
    xhatArray = [xhatArray xhat];
    % Simulate the system
    x = [1 dt;-k*dt/m 1-(b*dt/m)]*x+B*u;
    y = C*x;
    % Prediction
    F = [1 dt 0;
         -k*dt/m 1-(xhat(3)*dt/m) -dt*xhat(2)/m;
         0 0 1];
    xhat  = [xhat(1)+xhat(2)*dt;(-k*dt/m)*xhat(1)+(1-(xhat(3)*dt/m))*xhat(2)+dt*u/m;xhat(3)];
    Pmin  = F*Pplus*F' + Q;
    % Update
    K     = Pmin*H'*inv(H*Pmin*H' + R);
    Pplus = (eye(3)-K*H)*Pmin;
    xhat  = xhat + K*(y-H*xhat);
end

figure(1)
subplot(3,1,1)
plot(dt:dt:tf,xArray(1,:),'-b','LineWidth',3)
hold on;
plot(dt:dt:tf,xhatArray(1,:),':r','LineWidth',3)
legend('true position','estimated position')
set(gca,'FontSize',24)
grid on;
grid minor;
subplot(3,1,2)
plot(dt:dt:tf,xArray(2,:),'-b','LineWidth',3)
hold on;
plot(dt:dt:tf,xhatArray(2,:),':r','LineWidth',3)
legend('true velocity','estimated velocity')
set(gca,'FontSize',24)
grid on;
grid minor;
subplot(3,1,3)
plot(dt:dt:tf,[5*ones(1,5000) 6*ones(1,5000)],'-b','LineWidth',3)
hold on;
plot(dt:dt:tf,xhatArray(3,:),':r','LineWidth',3)
legend('true damping coefficient','estimated damping coefficient')
ylim([4.5 6.5])
set(gca,'FontSize',24)
grid on;
grid minor;