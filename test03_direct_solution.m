% Initialize camera object
num_cam = 3;
for i = 1:num_cam
    cam(i) = Camera();
end

% Set matrices for Projection matrix
K1 = cameraParams.IntrinsicMatrix; K2 = K1; K3 = K1;
R1 = cameraParams.RotationMatrices(:,:,1); R2 = cameraParams.RotationMatrices(:,:,2); R3 = cameraParams.RotationMatrices(:,:,3);
t1 = cameraParams.TranslationVectors(1,:)'; t2 = cameraParams.TranslationVectors(2,:)'; t3 = cameraParams.TranslationVectors(3,:)';

% Update Projection matrix and update G
cam(1).assign_KRt(K1, R1, t1);
cam(2).assign_KRt(K2, R2, t2);
cam(3).assign_KRt(K3, R3, t3);

for i = 1:num_cam
    cam(i).update_G();
end

% Update c_star
detection1 = [1900 2184-1200 400 1100];
detection2 = [2057 2184-758 500 2000];
detection3 = [1750 2184-830 500 1500];

cam(1).detection(detection1);
cam(2).detection(detection2);
cam(3).detection(detection3);
cam(1).update_c_star();
cam(2).update_c_star();
cam(3).update_c_star();

% Construct M for Optimization Problem
num_cam = 3;
M = zeros(6*num_cam, 10+num_cam);
for i = 1:num_cam
    M(6*(i-1)+1:6*i, 1:10) = cam(i).G;
    M(6*(i-1)+1:6*i, 10+i) = -cam(i).c_star;
end

% Solve w
[U, S, V] = svd(M);
S = 0.01*floor(100*S);
idx_min_singular = find(S, 1, 'last');
[idx_row, idx_col] = ind2sub(size(M), idx_min_singular);
w = V(:,idx_row);
v = w(1:10);
Q_adj = vech_inverse(v);
Q = adj_inverse(Q_adj);


function output = vech_inverse(v)
    dim = length(v);
    dim_output = floor(sqrt(dim*2));
    output = zeros(dim_output);
    count = 1;
    for i = 1:dim_output
        for j = i:dim_output
            output(j,i) = v(count);
            output(i,j) = v(count);
            count = count + 1;
        end
    end
end

function A = adj_inverse(A_adj)
    dim = length(A_adj);
    det_A = nthroot(det(A_adj),dim-1);
    A = -det_A * inv(A_adj);
end