%% Robotics I: Analysis - Control - Laboratory
% Semester 7 - Flow S

% Semester Project: Robotic Manipulator with 3 Rotational DOF
% Dimitris Dimos [031 17 165]


% Project Part B: Kinematic Simulation


%% function to simulate motion
% input: P_A (starting point), P_B (end point), T (motion duration) 
function simulate (P_A, P_B, T)

close all;  % close all currently opened figures

% translate input - trajectory end point coordinates
[x_A, y_A, z_A] = deal(P_A(1), P_A(2), P_A(3));  
[x_B, y_B, z_B] = deal(P_B(1), P_B(2), P_B(3));

if z_A ~= z_B
    disp('Error: Line Points must have same z-coordinate');
    return;
end

% robot links
l_0 = 0.0;
l_1 = 0.0;
l_2 = 4.0;
l_3 = 0.0;
l_4 = 14.0;
l_5 = 16.0;

% sampling period
dt = 0.001;  % 1kHz sampling frequency 

% trajectory polynomial coefficients
a0 = x_A;
a = [0, -3*(x_A - x_B)/T.^2, 2*(x_A - x_B)/T.^3];

% time axis definition
t = 0:dt:T;

% trajectory definition - end-effector position
x = a0 + a(1)*t + a(2)*t.^2 + a(3)*t.^3;
y = y_A + (x - x_A)*(y_B - y_A)/(x_B - x_A);
for i = 1:length(x)
    z(i) = z_A;
end

% trajectory definition - end-effector linear velocity
v_x = [0 diff(x)];
v_y = [0 diff(y)];
for i = 1:length(v_x)
    v_z(i) = 0;
end

%% PLOT: x - y linearity
figure(1)
xlabel('x(t)')
plot(x,y)
title('Linear Relationship y(t) - x(t)')
xlabel('x(t)')
ylabel('y = f(x(t))')
grid on

%% PLOT: desired position coordinates
figure(2)
subplot(3,1,1)
plot(t,x)
title('x - coordinate over time')
xlabel('time (sec)')
ylabel('x(t)')
grid on
subplot(3,1,2)
plot(t,y)
title('y - coordinate over time')
xlabel('time (sec)')
ylabel('y(t)')
grid on
subplot(3,1,3)
plot(t,z)
title('z - coordinate over time')
xlabel('time (sec)')
ylabel('z(t)')
grid on

%% PLOT: desired 3D trajectory
figure(3)
plot3(x,y,z)
xlim([x_A-1, x_B+1])
ylim([y_A-1, y_B+1])
zlim([z_A-1, z_B+1])
hold on
plot3(x_A,y_A,z_A,'O')
plot3(x_B,y_B,z_B,'O')
xlabel('x(t)')
ylabel('y(t)')
zlabel('z(t)')
grid on

%% PLOT: desired velocity components
figure(4)
subplot(3,1,1)
plot(t,v_x)
title('v_x - component over time')
xlabel('time (sec)')
ylabel('v_x(t)')
grid on
subplot(3,1,2)
plot(t,v_y)
title('v_y - component over time')
xlabel('time (sec)')
ylabel('v_y(t)')
grid on
subplot(3,1,3)
plot(t,v_z)
title('v_z - component over time')
xlabel('time (sec)')
ylabel('v_z(t)')
grid on

%% INVERSE KINEMATIC MODEL
nom_3 = (x(:) - l_1).^2 + y(:).^2 + (z(:) + l_0).^2 - l_2.^2 - l_4.^2 - l_5.^2;
den_3 = 2 * l_4 * l_5;
q_3 = acos(nom_3 / den_3);

nom_2 = y(:).*l_5.*sin(q_3) + (l_4 + l_5.*cos(q_3)).*sqrt( (l_4 + l_5.*cos(q_3)).^2 + (l_5.*sin(q_3) ).^2 - y(:).^2);
den_2 = (l_5 .* sin(q_3)) .^ 2 + (l_4 + l_5 * cos(q_3)).^2;
q_2 = acos(nom_2./den_2);

nom_1 = (l_4.*cos(q_2) + l_5.*cos(q_2 + q_3)).*z(:) + l_2.*sqrt((l_4.*cos(q_2) + l_5.*cos(q_2+q_3)).^2 + (l_1 + l_2).^2 - z(:).^2);
den_1 = (l_4.*cos(q_2) + l_5.*cos(q_2 + q_3)).^2 + l_2.^2;
q_1 = acos(nom_1 ./ den_1);

%% PLOT: desired joint angles
figure(5)
subplot(3,1,1)
plot(t,q_3(:))
title('q_3 - angle over time')
xlabel('time (sec)')
ylabel('q_3(t)')
grid on
subplot(3,1,2)
plot(t,q_2(:))
title('q_2 - angle over time')
xlabel('time (sec)')
ylabel('q_2(t)')
grid on
subplot(3,1,3)
plot(t,q_1(:))
title('q_1 - angle over time')
xlabel('time (sec)')
ylabel('q_1(t)')
grid on

%% INVERSE DIFFERENTIAL MODEL

w_1 = [0; diff(q_1)];
w_2 = [0; diff(q_2)];
w_3 = [0; diff(q_3)];

%% PLOT: desired joint angle velocity
figure(6)
subplot(3,1,1)
plot(t,w_3(:))
title('w_3 - angle velocity over time')
xlabel('time (sec)')
ylabel('w_3(t)')
grid on
subplot(3,1,2)
plot(t,w_2(:))
title('w_2 - angle velocity over time')
xlabel('time (sec)')
ylabel('w_2(t)')
grid on
subplot(3,1,3)
plot(t,w_1(:))
title('w_1 - angle velocity over time')
xlabel('time (sec)')
ylabel('w_1(t)')
grid on

%% FORWARD KINEMATIC MODEL

points = length(t);

% joint 1 position
j1_x(1:points) = 0;
j1_y(1:points) = 0;
j1_z(1:points) = 0;

% joint 2 position
j2_x = l_2 .* cos(q_1);
j2_y = l_2 .* sin(q_1);
j2_z(1:points) = 0;

% joint 3 position
j3_x = l_1 + cos(q_1) .* l_2 + l_3 .* sin(q_1) + l_4 .* cos(q_2) .* sin(q_1);
j3_y = l_4 .* sin(q_2);
j3_z = -l_0 + l_2 .* sin(q_1) - l_3 .* cos(q_1) - l_4 .* cos(q_1) .* cos(q_2);



%% KINEMATIC SIMULATION
figure(7)

plot3(x,y,z, 'gs'); % trajectory
hold on
dt_a = 1000; 
plot([0], [0], 'o')    % axis origin

%% animation
for i = 1:dt_a:points
    pause(0.1);     % pause motion
    
    % joint 1
    plot3([j1_x(i)],[j1_y(i)],[j1_z(i)],'bs');
    
    % link from joint 1 to joint 2
    plot3([j1_x(i),j2_x(i)],[j1_y(i),j2_y(i)],[j1_z(i),j2_z(i)], 'k');
    
    % joint 2
    plot3([j2_x(i)],[j2_y(i)],[j2_z(i)],'bs');
    
    % link from joint 2 to joint 3
    plot3([j2_x(i),j3_x(i)],[j2_y(i),j3_y(i)],[j2_z(i),j3_z(i)], 'k');
    
    % joint 3
    plot3([j3_x(i)],[j3_y(i)],[j3_z(i)],'bs');
    
    % link from joint 3 to end-effector
    plot3([j3_x(i),x(i)],[j3_y(i),y(i)],[j3_z(i),z(i)], 'k');
    
    % end-effector
    plot3([x(i)],[y(i)],[z(i)],'bx');
end

xlabel('x - axis');
ylabel('y - axis');
zlabel('z - axis');
grid on

end
