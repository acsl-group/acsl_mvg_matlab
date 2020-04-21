% Initialize camera object
cam1 = Camera();

% Set matrices for Projection matrix
K = eye(3);
R = diag([0.5,0.5,1]);
t = [0 0 0]';

% Update Projection matrix and update G
cam1.assign_KRt(K, R, t);
cam1.update_G();

% Update c_star
xc = 0.3; yc = 0.2; w = 0.8; h = 0.5;
cam1.detection([xc, yc, w, h]);
cam1.update_c_star();

% Result
cam1.G
cam1.c_star