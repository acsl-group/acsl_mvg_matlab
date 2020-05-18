% Projection matrix
K = sym('K%d%d', [3 3]);
R = sym('R%d%d', [3 3]);
t = sym('t%d', [3 1]);
P = K * [R t];
% P = sym('P%d%d', [3 4]);

% C in Conic equation (2D ellipse)
syms xc yc w h
assume([xc, yc, w, h],'positive');
C = [4*w^2 0 -4*w^2*xc ;
    0 4*h^2 -4*h^2*yc ;
    -4*w^2*xc -4*h^2*yc 4*w^2*xc^2+4*h^2*yc^2-w^2*h^2];
C_adj = adjoint(C);
C_adj = -C_adj / C_adj(3,3);
c_star = vech(C_adj);
adj_inverse(C_adj)

% D and E (Eq. (4))
D = zeros(6,9);
D(1,1) = 1; D(2,2) = 1; D(3,3) = 1; D(4,5) = 1; D(5,6) = 1; D(6,9) = 1;
E = zeros(16,10);
E(1,1) = 1; E(2,2) = 1; E(3,3) = 1; E(4,4) = 1; E(5,2) = 1; E(6,5) = 1;
E(7,6) = 1; E(8,7) = 1; E(9,3) = 1; E(10,6) = 1; E(11,8) = 1; E(12,9) = 1;
E(13,4) = 1; E(14,7) = 1; E(15,9) = 1; E(16,10) = 1;

% G (Eq. (4))
G = D * kron(P,P) * E;



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
    A = det_A * inv(A_adj) / 1i;
end