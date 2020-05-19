%% Ground truth Q
Q = [36 0 0 -36 ; 0 9 0 -18 ; 0 0 4 -20 ; -36 -18 -20 136];     % (x-1)^2 + (y-2)^2 / 4 + (z-5)^2 / 9 = 1

% Obtain rotation transform R
A = Q(1:3,1:3);
[V_eig, d_eig] = eig(A, 'vector');      % A = V_eig * diag(d_eig) * V_eig'
R = inv(V_eig);                         % A = R' * diag(d_eig) * R

% Obtain translation
trans = (R' * diag(d_eig)) \ Q(1:3,4);

% Normalization
normalize_factor = trans' * diag(d_eig) * trans - Q(4,4);
d_eig = d_eig / normalize_factor;

% Axis length
length_x = nthroot(d_eig(1), -2);
length_y = nthroot(d_eig(2), -2);
length_z = nthroot(d_eig(3), -2);

[x_transformed,y_transformed,z_transformed] = ellipsoid(0,0,0,length_x,length_y,length_z,30);
x = zeros(size(x_transformed));
y = zeros(size(y_transformed));
z = zeros(size(z_transformed));

num_plot = length(x_transformed);
for i = 1:num_plot
    for j = 1:num_plot
        T = [R trans; zeros(1,3) 1];
        coord = T \ [x_transformed(i,j) y_transformed(i,j) z_transformed(i,j) 1]';
        coord = coord / coord(4);  % normalize
        x(i,j) = coord(1);
        y(i,j) = coord(2);
        z(i,j) = coord(3);
    end
end

%% Camera parameters
load('cameraParams.mat')

% Load intrinsic parameter matrix
K = cameraParams.IntrinsicMatrix';

% Initialize camera object
num_cam = 3;
for i = 1:num_cam
    cam(i) = Camera();
end

% Euler angles (yaw, pitch, roll)
% Drone - Camera
eul_DC = [-pi/2 -pi/2 pi/2]';

% Euler angles (yaw, pitch, roll) and Position of Drones
% World - Drone
eul1 = [-pi/4, pi - acos(sqrt(2/3)), 0]';
eul2 = [pi/4, pi - acos(sqrt(2/3)), 0]';
eul3 = [0 0 0]';
position1 = [20 -20 20]';
position2 = [20 20 20]';
position3 = [100 0 0]';

% Compute Extrinsic matrix parameters
R_DC = Rotation_matrix(eul_DC(1), eul_DC(2), eul_DC(3));
Rc1 = Rotation_matrix(eul1(1), eul1(2), eul1(3));
Rc2 = Rotation_matrix(eul2(1), eul2(2), eul2(3));
Rc3 = Rotation_matrix(eul3(1), eul3(2), eul3(3));
Transform1 = [Rc1 position1 ; zeros(1,3) 1] * [R_DC zeros(3,1) ; zeros(1,3) 1];
Transform2 = [Rc2 position2 ; zeros(1,3) 1] * [R_DC zeros(3,1) ; zeros(1,3) 1];
Transform3 = [Rc3 position3 ; zeros(1,3) 1] * [R_DC zeros(3,1) ; zeros(1,3) 1];
R1 = Transform1(1:3,1:3)';
R2 = Transform2(1:3,1:3)';
R3 = Transform3(1:3,1:3)';
t1 = -R1 * Transform1(1:3,4);
t2 = -R2 * Transform2(1:3,4);
t3 = -R3 * Transform3(1:3,4);

% Update Projection matrix and update G
cam(1).assign_KRt(K, R1, t1);
cam(2).assign_KRt(K, R2, t2);
cam(3).assign_KRt(K, R3, t3);

for i = 1:num_cam
    cam(i).update_G();
end

%% Detection (Obtain xc, yc, w, h)
% Camera 1
x_cam1 = zeros(num_plot, num_plot);
y_cam1 = zeros(num_plot, num_plot);
for i = 1:num_plot
    for j = 1:num_plot
        result_temp = K*[R1 t1] * [x(i,j) y(i,j) z(i,j) 1]';
        result_temp = result_temp / result_temp(3);
        x_cam1(i,j) = result_temp(1);
        y_cam1(i,j) = result_temp(2);
    end
end

% Camera 2
x_cam2 = zeros(num_plot, num_plot);
y_cam2 = zeros(num_plot, num_plot);
for i = 1:num_plot
    for j = 1:num_plot
        result_temp = K*[R2 t2] * [x(i,j) y(i,j) z(i,j) 1]';
        result_temp = result_temp / result_temp(3);
        x_cam2(i,j) = result_temp(1);
        y_cam2(i,j) = result_temp(2);
    end
end

% Camera 3
x_cam3 = zeros(num_plot, num_plot);
y_cam3 = zeros(num_plot, num_plot);
for i = 1:num_plot
    for j = 1:num_plot
        result_temp = K*[R3 t3] * [x(i,j) y(i,j) z(i,j) 1]';
        result_temp = result_temp / result_temp(3);
        x_cam3(i,j) = result_temp(1);
        y_cam3(i,j) = result_temp(2);
    end
end

% Plot point detection results
figure(1); plot(reshape(x_cam1,[],1), reshape(y_cam1,[],1),'.');
xlim([1 4608]); ylim([1 2184]); ax = gca; ax.YDir = 'reverse';
title('Camera 1')
figure(2); plot(reshape(x_cam2,[],1), reshape(y_cam2,[],1),'.');
xlim([1 4608]); ylim([1 2184]); ax = gca; ax.YDir = 'reverse';
title('Camera 2')
figure(3); plot(reshape(x_cam3,[],1), reshape(y_cam3,[],1),'.');
xlim([1 4608]); ylim([1 2184]); ax = gca; ax.YDir = 'reverse';
title('Camera 3')

% Consider bounding box of each image
x_min_cam1 = min(x_cam1,[],'all'); y_min_cam1 = min(y_cam1,[],'all');
x_max_cam1 = max(x_cam1,[],'all'); y_max_cam1 = max(y_cam1,[],'all');

x_min_cam2 = min(x_cam2,[],'all'); y_min_cam2 = min(y_cam2,[],'all');
x_max_cam2 = max(x_cam2,[],'all'); y_max_cam2 = max(y_cam2,[],'all');

x_min_cam3 = min(x_cam3,[],'all'); y_min_cam3 = min(y_cam3,[],'all');
x_max_cam3 = max(x_cam3,[],'all'); y_max_cam3 = max(y_cam3,[],'all');

% xc, yc, w, h
xc_cam1 = (x_min_cam1 + x_max_cam1) / 2;
yc_cam1 = (y_min_cam1 + y_max_cam1) / 2;
w_cam1 = y_max_cam1 - y_min_cam1;
h_cam1 = x_max_cam1 - x_min_cam1;

xc_cam2 = (x_min_cam2 + x_max_cam2) / 2;
yc_cam2 = (y_min_cam2 + y_max_cam2) / 2;
w_cam2 = y_max_cam2 - y_min_cam2;
h_cam2 = x_max_cam2 - x_min_cam2;

xc_cam3 = (x_min_cam3 + x_max_cam3) / 2;
yc_cam3 = (y_min_cam3 + y_max_cam3) / 2;
w_cam3 = y_max_cam3 - y_min_cam3;
h_cam3 = x_max_cam3 - x_min_cam3;


function R = Rotation_matrix(yaw, pitch, roll)
    Rz = [cos(yaw) -sin(yaw) 0 ; sin(yaw) cos(yaw) 0 ; 0 0 1];
    Ry = [cos(pitch) 0 sin(pitch) ; 0 1 0 ; -sin(pitch) 0 cos(pitch)];
    Rx = [1 0 0 ; 0 cos(roll) -sin(roll) ; 0 sin(roll) cos(roll)];
    R = Rz * Ry * Rx;
end