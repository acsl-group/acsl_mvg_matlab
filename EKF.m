% 2020-07-30
classdef EKF < handle
    properties
        % functions
        f, h, Jacobian_F, Jacobian_H;
        % covariance matrices
        P, Q, R;
        % state
        x_pre;
    end
    methods
        %% function area
        function obj = EKF(P, Q, R, f, Jacobian_F, h, Jacobian_H, state_init)
            % matrix init
            obj.P = P;
            obj.Q = Q;
            obj.R = R;
            % function init
            obj.f = f;
            obj.h = h;
            obj.Jacobian_F = Jacobian_F;
            obj.Jacobian_H = Jacobian_H;
            
            obj.x_pre = state_init;
        end
        
        function x_hat = Estimate(obj, u, z)
            x_size = size(obj.x_pre, 1);
            arguments = num2cell([obj.x_pre' u']);
            F = obj.Jacobian_F(arguments{:});
            H = obj.Jacobian_H(arguments{:});
            % Prediction
            state_hat_temp = obj.f(arguments{:});
            obj.P = F * obj.P * F' + obj.Q;
            K = obj.P * H' / (H*obj.P*H' + obj.R);
            % Innovation
            Inno = z - obj.h(arguments{:});
            % Correction
            x_hat = state_hat_temp + K * Inno;
            obj.P = (eye(x_size) - K*H) * obj.P * (eye(x_size) - K*H)' + K*obj.R*K';
            
            obj.x_pre = x_hat;
        end
    end
end
