clear all
load('test09_data.mat')

%% System model
syms state [20 1] real
syms dt real
% State equation
f(state) = [eye(10) dt*eye(10); zeros(10) eye(10)] * state;

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

state_init = zeros(20, 1);
Est1 = EKF(eye(20), eye(20), eye(18), f, Jacobian_F, h, Jacobian_H, state_init);


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