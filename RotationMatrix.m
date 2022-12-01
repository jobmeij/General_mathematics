%% Rotation matrix examples to visualize what happens
clear all; close all; clc

%% Starting in two dimensions

vect2d = [1; 0]

% Create 2D rotation matrix
t = -45     % Angle of rotation 
rot2d = [cosd(t) sind(t); -sind(t) cosd(t)]
% Matrix above corresponds to the complex number like
cosd(t) + i*sind(t) 
exp(i*t*(pi/180))

vect2d_rot = rot2d*vect2d

% Figure of vector before and after rotation matrix
figure()
hold on
plot([0, vect2d(1)],[0, vect2d(2)])
plot(vect2d(1),vect2d(2),'*')
%
plot([0, vect2d_rot(1)],[0, vect2d_rot(2)])
plot(vect2d_rot(1),vect2d_rot(2),'*')
hold off
grid on
xlabel('x')
ylabel('y')

%% three dimension rotation matrices 
theta_x = 90
theta_y = 90
theta_z = 90

vect3d = [1; 0; 0]

Rx =    [1 0 0;
         0 cosd(theta_x) -sind(theta_x);
         0 sind(theta_x) cosd(theta_x)]

Ry =    [cosd(theta_y) 0 sind(theta_y);
         0 1 0;
         -sind(theta_y) 0 cosd(theta_y)]

Rz =    [cosd(theta_z) -sind(theta_z) 0;
         sind(theta_z) cosd(theta_z) 0;
         0 0 1]

% Compute rotated vectors
vect3d_rot_x = Rx*vect3d
vect3d_rot_y = Ry*vect3d
vect3d_rot_z = Rz*vect3d

% Figure of vector before and after rotation matrix
figure()
hold on
plot3([0, vect3d(1)],[0, vect3d(2)],[0, vect3d(3)])
plot3(vect3d(1),vect3d(2),vect3d(3),'*')
%
plot3([0, vect3d_rot_x(1)],[0, vect3d_rot_x(2)],[0, vect3d_rot_x(3)])
plot3(vect3d_rot_x(1),vect3d_rot_x(2),vect3d_rot_x(3),'*')
hold off
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('X-axis rotation matrix')

% Figure of vector before and after rotation matrix
figure()
hold on
plot3([0, vect3d(1)],[0, vect3d(2)],[0, vect3d(3)])
plot3(vect3d(1),vect3d(2),vect3d(3),'*')
%
plot3([0, vect3d_rot_y(1)],[0, vect3d_rot_y(2)],[0, vect3d_rot_y(3)])
plot3(vect3d_rot_y(1),vect3d_rot_y(2),vect3d_rot_y(3),'*')
hold off
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Y-axis rotation matrix')

% Figure of vector before and after rotation matrix
figure()
hold on
plot3([0, vect3d(1)],[0, vect3d(2)],[0, vect3d(3)])
plot3(vect3d(1),vect3d(2),vect3d(3),'*')
%
plot3([0, vect3d_rot_z(1)],[0, vect3d_rot_z(2)],[0, vect3d_rot_z(3)])
plot3(vect3d_rot_z(1),vect3d_rot_z(2),vect3d_rot_z(3),'*')
hold off
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Z-axis rotation matrix')


%% General rotations in 3D
syms sin_psi sin_gamma sin_phi cos_psi cos_gamma cos_phi
R_EB_x =   [cos_psi, 0, sin_psi;
            0, 1, 0;
            -sin_psi, 0, cos_psi]

R_EB_y =   [1, 0, 0;
            0, cos_gamma, -sin_gamma
            0, sin_gamma, cos_gamma]

R_EB_z =   [cos_phi, -sin_phi, 0;
            sin_phi, cos_phi, 0;
            0, 0, 1]

% Generating main rotation matrix
R_EB = R_EB_x*R_EB_y*R_EB_z







