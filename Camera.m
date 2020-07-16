classdef Camera < handle
    % 카메라 정보를 받아와서 Optimization problem을 풀기 위한 c_star, G 준비
    % P: Projection matrix
    % xc, yc, w, h => C => c_star
    % P, D, E => G
    
    properties
        K, R, t, P, D, E, G, H, T
        xc, yc, w, h
        C, c_star, C2, c2_star
        C_func, c_star_func, C2_func, c2_star_func
    end
    
    methods
        function obj = Camera()
            % Projection matrix
            obj.K = eye(3); obj.R = eye(3); obj.t = [0 0 0]';
            obj.P = obj.K*[obj.R obj.t];
            
            obj.D = zeros(6,9);
            obj.D(1,1) = 1; obj.D(2,2) = 1; obj.D(3,3) = 1; obj.D(4,5) = 1;
            obj.D(5,6) = 1; obj.D(6,9) = 1;
            
            obj.E = zeros(16,10);
            obj.E(1,1) = 1; obj.E(2,2) = 1; obj.E(3,3) = 1; obj.E(4,4) = 1;
            obj.E(5,2) = 1; obj.E(6,5) = 1; obj.E(7,6) = 1; obj.E(8,7) = 1;
            obj.E(9,3) = 1; obj.E(10,6) = 1; obj.E(11,8) = 1;
            obj.E(12,9) = 1; obj.E(13,4) = 1; obj.E(14,7) = 1;
            obj.E(15,9) = 1; obj.E(16,10) = 1;
            
            obj.G = obj.D * kron(obj.P, obj.P) * obj.E;
            obj.xc = []; obj.yc = []; obj.w = []; obj.h = [];
            
            % Define functions for generate C and c_star
            syms xc yc w h
            assume([xc, yc, w, h],'positive');
            C = [4*h^2 0 -4*h^2*xc ;
                0 4*w^2 -4*w^2*yc ;
                -4*h^2*xc -4*w^2*yc 4*h^2*xc^2+4*w^2*yc^2-w^2*h^2];
            C2 = [4/w^2 0 0 ; 
                0 4/h^2 0 ;
                0 0 -1];
            C_adj = adjoint(C);
            C2_adj = adjoint(C2);
            C_adj = -C_adj / C_adj(3,3);
            C2_adj = -C2_adj / C2_adj(3,3);
            c_star = vech(C_adj);
            c2_star = vech(C2_adj);
            
            C_func = symfun(C, [xc, yc, w, h]);
            C2_func = symfun(C2, [xc, yc, w, h]);
            obj.C_func = matlabFunction(C_func);
            obj.C2_func = matlabFunction(C2_func);
            c_star = symfun(c_star, [xc, yc, w, h]);
            c2_star = symfun(c2_star, [xc, yc, w, h]);
            obj.c_star_func = matlabFunction(c_star);
            obj.c2_star_func = matlabFunction(c2_star);
        end
        
        function assign_P(obj, P)
            obj.P = P;
        end
        
        function assign_KRt(obj, K, R, t)
            obj.K = K; obj.R = R; obj.t = t;
            obj.P = obj.K*[obj.R obj.t];
        end
        
        function update_G(obj)
            obj.G = obj.D * kron(obj.P, obj.P) * obj.E;
        end
        
        function detection(obj, rectangle)
            % rectangle: [xc, yc, w, h]'
            obj.xc = rectangle(1);
            obj.yc = rectangle(2);
            obj.w = rectangle(3);
            obj.h = rectangle(4);
        end
        
        function update_C(obj)
            obj.C = obj.C_func(obj.xc, obj.yc, obj.w, obj.h);
        end
        
        function update_c_star(obj)
            obj.c_star = obj.c_star_func(obj.xc, obj.yc, obj.w, obj.h);
        end
        
        function update_c2_star(obj)
            obj.c2_star = obj.c2_star_func(obj.xc, obj.yc, obj.w, obj.h);
        end
        
        function update_H(obj)
            h_ellipse = sqrt((obj.w * 0.5)^2 + (obj.h * 0.5)^2);
            obj.H = zeros(3);
            obj.H(1:2, 1:2) = h_ellipse * eye(2);
            obj.H(1:2, 3) = [obj.xc obj.yc]';
            obj.H(3, 3) = 1;
        end
        
        function update_T(obj, Q_adj)
            obj.T = eye(4);
            obj.T(1:3, 4) = Q_adj(1:3, 4);
        end
    end
end

function output = vech(A)
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