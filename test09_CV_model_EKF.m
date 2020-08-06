clear all
load('test09_data.mat')
load('bbox.mat')

%% System model
syms state [20 1] real

% State equation
f(state) = [eye(10) eye(10); zeros(10) eye(10)] * state;

% Measurement equation
Q_star = vech_inverse_sym([eye(10) zeros(10)] * state);
C_star1 = P1 * Q_star * P1';
C_star2 = P2 * Q_star * P2';
C_star3 = P3 * Q_star * P3';
c_star1 = vech_sym(C_star1);
c_star2 = vech_sym(C_star2);
c_star3 = vech_sym(C_star3);
h(state) = simplify([c_star1' c_star2' c_star3']');

%% Jacobians
Jacobian_F = jacobian(f, state');
Jacobian_H = jacobian(h, state');

%% Define EKF
f = matlabFunction(f);
Jacobian_F = matlabFunction(Jacobian_F);
h = matlabFunction(h);
Jacobian_H = matlabFunction(Jacobian_H);

x_hat = zeros(20, 1);
Est1 = EKF(eye(20), eye(20), eye(18), f, Jacobian_F, h, Jacobian_H, x_hat);

%% Define Cameras
for i = 1:3
    cam(i) = Camera();
end
cam(1).assign_P(P1);
cam(2).assign_P(P2);
cam(3).assign_P(P3);

%% Simulation
num_iteration = size(bbox.cam1, 2);
num_plot = 30;
figure();
for i = 1:num_iteration
    cam(1).detection(bbox.cam1(:,i));
    cam(2).detection(bbox.cam2(:,i));
    cam(3).detection(bbox.cam3(:,i));
    
    cam(1).update_c_star();
    cam(2).update_c_star();
    cam(3).update_c_star();
    
    msr = [cam(1).c_star' cam(2).c_star' cam(3).c_star']';
    ctrl = [];
    
    x_hat = Est1.Estimate(ctrl, msr);
    Q_star = vech_inverse(x_hat(1:10));
    try
        [x,y,z] = ellipsoid_from_Q(Q_star, num_plot);
        surf(x,y,z); hold on;
        drawnow;
    catch
        
    end
end




function output = vech_sym(A)
% Serializes the elements of the lower triangular part
    dim = length(A);
    dim_output = dim*(dim+1)/2;
    output = sym(zeros(dim_output, 1));
    count = 1;
    for i = 1:dim
        for j = i:dim       % from "i" to dim
            output(count) = A(j,i);
            count = count + 1;
        end
    end
end

function output = vech_inverse_sym(v)
    dim = length(v);
    dim_output = floor(sqrt(dim*2));
    output = sym(zeros(dim_output));
    count = 1;
    for i = 1:dim_output
        for j = i:dim_output
            output(j,i) = v(count);
            output(i,j) = v(count);
            count = count + 1;
        end
    end
end

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

function [x, y, z] = ellipsoid_from_Q(Q, num_plot)
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

    [x_transformed,y_transformed,z_transformed] = ellipsoid(0,0,0,length_x,length_y,length_z,num_plot);
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

end