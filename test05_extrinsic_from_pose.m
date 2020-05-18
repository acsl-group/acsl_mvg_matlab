% Building the Extrinsic Matrix from Camera Pose
% Link: http://ksimek.github.io/2012/08/22/extrinsic/
% Yaw, pitch, and roll to a rotation matrix
% Link: http://planning.cs.uiuc.edu/node102.html

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
eul_DC = [pi/2 0 pi/2]';

% Euler angles (yaw, pitch, roll) and Position of Drones
% World - Drone
eul1 = [-pi/4, pi - acos(sqrt(2/3)), 0]';
eul2 = [pi/4, pi - acos(sqrt(2/3)), 0]';
eul3 = [0 0 0]';
position1 = [10 -10 10]';
position2 = [10 10 10]';
position3 = [10 0 0]';

% Compute Extrinsic matrix parameters
R_DC = Rotation_matrix(eul_DC(1), eul_DC(2), eul_DC(3));
Rc1 = Rotation_matrix(eul1(1), eul1(2), eul1(3));
Rc2 = Rotation_matrix(eul2(1), eul2(2), eul2(3));
Rc3 = Rotation_matrix(eul3(1), eul3(2), eul3(3));
Transform1 = [R_DC zeros(3,1) ; zeros(1,3) 1] * [Rc1 position1 ; zeros(1,3) 1];
Transform2 = [R_DC zeros(3,1) ; zeros(1,3) 1] * [Rc2 position2 ; zeros(1,3) 1];
Transform3 = [R_DC zeros(3,1) ; zeros(1,3) 1] * [Rc3 position3 ; zeros(1,3) 1];
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


function R = Rotation_matrix(yaw, pitch, roll)
    Rz = [cos(yaw) -sin(yaw) 0 ; sin(yaw) cos(yaw) 0 ; 0 0 1];
    Ry = [cos(pitch) 0 sin(pitch) ; 0 1 0 ; -sin(pitch) 0 cos(pitch)];
    Rx = [1 0 0 ; 0 cos(roll) -sin(roll) ; 0 sin(roll) cos(roll)];
    R = Rz * Ry * Rx;
end