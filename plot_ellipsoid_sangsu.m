% Run this code after running test03 code

% Assume there is transformation T such that x = T * x_transformed
[V_eig, d_eig] = eig(Q, 'vector');                  % Q = V_eig * diag(d_eig) * V_eig'
[d_eig, ind] = sort(d_eig, 'descend');
V_eig = V_eig(:, ind);
normalize_factor = -d_eig(4);
Q_transformed = diag(d_eig') / normalize_factor;
T = inv(V_eig' * sqrt(normalize_factor));           % Q = inv(T)' * Q_transformed * inv(T)

length_x = nthroot(Q_transformed(1,1), -2);
length_y = nthroot(Q_transformed(2,2), -2);
length_z = nthroot(Q_transformed(3,3), -2);

[x_transformed,y_transformed,z_transformed] = ellipsoid(0,0,0,length_x,length_y,length_z,30);
x = zeros(size(x_transformed));
y = zeros(size(y_transformed));
z = zeros(size(z_transformed));

num_plot = length(x_transformed);
for i = 1:num_plot
    for j = 1:num_plot
        coord = T * [x_transformed(i,j) y_transformed(i,j) z_transformed(i,j) 1]';
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