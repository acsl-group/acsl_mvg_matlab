% Q = diag([4,1,9,-1]);
% Transform1 = [cos(pi/4) -sin(pi/4) 0 0 ; sin(pi/4) cos(pi/4) 0 0 ; 0 0 1 0 ; 0 0 0 1];
% Transform2 = [cos(pi/6) 0 -sin(pi/6) 0 ; 0 1 0 0 ; sin(pi/6) 0 cos(pi/6) 0 ; 0 0 0 1];
% Transform3 = [1 0 0 3 ; 0 1 0 7 ; 0 0 1 5 ; 0 0 0 1];
% Q = inv(Transform1)' * Q * inv(Transform1);
% Q = inv(Transform2)' * Q * inv(Transform2);
% Q = inv(Transform3)' * Q * inv(Transform3);

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

figure
surf(x,y,z)
axis equal
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');