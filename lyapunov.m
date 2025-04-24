function [final_EX2]=mc(a,b,g,r,T,Nh,Nk,N,x0)
% Approximates the functional Phi(x) using a Monte Carlo method based on the fully discrete SPDE approximation.
    h=1/(Nh+1);
    k=T/Nk;

    mini_chol = cholesky_matrix_h_interval();
    A = a*stiffness_matrix(h);
    M = mass_matrix(h);
    level_EX2=0;
    parfor m = 1:N
        level_EX2 = level_EX2 + sample_HSHE_end_norm(A, M, T, k, h,b,x0,mini_chol,g,r);
    end
    final_EX2 = level_EX2/N;
end

function [x_val] = sample_HSHE_end_norm(A, M, T, tau, h,b,x0,mini_chol,g,r)
% Generates one sample path of the fully discrete SPDE approximation and computes the corresponding sample value for the functional Phi(x).
% Assumes that the noise is white in space and spatial dimension d = 1.
    nodes = h^(-1)-1;
    xgrid = linspace(0,1,nodes+2);
    timesteps = T/tau;
    x = x0(xgrid);
    x = x(2:end-1)';
    rr = x'*M*x*r^2*tau;
    for j = 1:timesteps
        W = generate_wiener_matrix(mini_chol,h,tau);
        x = (M + tau*A) \ (M*x + b*W*x);
        rr = rr + x'*M*x*r^2*tau;
    end
    x_val = rr + g^2*x'*M*x - x'*M*x*r^2*tau;
end

function [W] = generate_wiener_matrix(M_chol,h,tau)
% Generates a matrix representation of the spatial projection of the Wiener process increment used in the SPDE discretization.
    W = sparse(1/h-1,1/h-1);
    R = sqrt(h*tau/60)*M_chol*randn(3,1/h-2);
    W = spdiags([R(1,:) 0]',0,W);
    W = spdiags([[R(2,:)'; 0] [0; R(2,:)']],[-1,1],W);
    W = W + spdiags([0 R(3,:)]',0,1/h-1,1/h-1);
    W(1,1) = W(1,1) + sqrt(h*tau/60)*M_chol(1,1)*randn;
    W(1/h-1,1/h-1) = W(1/h-1,1/h-1) + sqrt(h*tau/60)*M_chol(1,1)*randn;
end

function A = stiffness_matrix(h)
% Computes the stiffness matrix A_h corresponding to the finite element discretization of the operator A = -Delta.
    A = h^(-1)*(2*spdiags(ones(h^(-1)-1,1),0,h^(-1)-1,h^(-1)-1) ...
        - spdiags(ones(h^(-1)-1,1),1,h^(-1)-1,h^(-1)-1) ...
            - spdiags(ones(h^(-1)-1,1),-1,h^(-1)-1,h^(-1)-1));
end

function M = mass_matrix(h)
% Computes the mass matrix M_h corresponding to the finite element discretization.
    M = h/6*(4*spdiags(ones(h^(-1)-1,1),0,h^(-1)-1,h^(-1)-1) ...
        + spdiags(ones(h^(-1)-1,1),1,h^(-1)-1,h^(-1)-1) ...
        + spdiags(ones(h^(-1)-2,1),-1,h^(-1)-1,h^(-1)-1));
end

function [M_chol] = cholesky_matrix_h_interval()
% Computes a Cholesky factor used in generate_wiener_matrix.
    A = [12 3 2; 3 2 3; 2 3 12];
    M_chol = chol(A);
end

function L = compute_lyapunov(a,b,g,r,T,Nh,Nk)
% Computes the matrix representation L^N_{h,tau} of the fully discrete approximation L^N_{h,tau} of the Lyapunov equation solution L(T).
    h=1/(Nh+1);
    k=T/Nk;

    S = a*stiffness_matrix(h);
    M = mass_matrix(h);

    L=g^2*M\spdiags(ones(Nh,1),0,Nh,Nh);
    B_id = sub2ind([Nh,Nh],2:Nh-1,2:Nh-1);
    B_iu = sub2ind([Nh,Nh],1:Nh-1,2:Nh);
    B_il = sub2ind([Nh,Nh],2:Nh,1:Nh-1);
    B = sparse(Nh,Nh);
    for i=1:Nk
        if Nh == 1
            B(1,1) = 2*h/5;
        else
            L_md = diag(L); L_ud = diag(L,1); L_ld = diag(L,-1);
            B(B_id) = h*(L_md(2:Nh-1)'*2/5+L_ld(2:Nh-1)'/10+L_ud(1:Nh-2)'/10+ ...
                     L_md(1:Nh-2)'/30+L_md(3:Nh)'/30);
            B(1,1) = h*(L(1,1)*2/5+L(1,2)/10+L(2,2)/30);
            B(Nh,Nh) = h*(L(Nh,Nh)*2/5+L(Nh,Nh-1)/10+L(Nh-1,Nh-1)/30);
            B(B_iu) = h*(L_md(1:Nh-1)'/20+L_md(2:Nh)'/20+L_ud(1:Nh-1)'/15);
            B(B_il) = B(B_iu);
        end
        RHS=M*L*M+k*(r^2*M+b^2*B);
        L = mldivide((M + k*S),RHS);
        L = mldivide((M + k*S),L');
    end
end

function Phi = compute_lyapunov_functional(L,u0)
% Computes the approximation Phi^L_{h,tau}(x) = <L^N_{h,tau} I_h x, I_h x> using the output of compute_lyapunov.
    Nh=length(L);
    h=1/(Nh+1);

    M=mass_matrix(h);

    X00=u0(linspace(h,1-h,Nh))';
    c=M*X00;
    Phi = 0;
    for i=1:Nh
        for j=1:Nh
            Phi = Phi + L(i,j)*c(i)*c(j);
        end
    end
end

% -------------------------------------------------------------------------
% Demo Usage
% -------------------------------------------------------------------------

% Define parameters for the SPDE and functional
a = 0.05; % Coefficient for the Laplacian term (-a*Delta)
b = 0.65; % Coefficient for the multiplicative noise term (b*X*dW)
g = 1.0;  % Coefficient for the terminal condition term (||g*X(T)||^2)
r = 0.0;  % Coefficient for the integrated term (int ||r*X(t)||^2 dt)
T = 1.0;  % Final time

% Define discretization parameters
Nh = 31;  % Number of interior spatial grid points (Nh+1 intervals)
Nk = Nh^2; % Number of time steps (Nk+1 points), chosen s.t. tau = h^2

% Define Monte Carlo parameters
N_mc = 1024; % Number of Monte Carlo samples

% Define the initial condition function handle
x0 = @(x) x.*(x < 0.5) + (1-x).*(x >= 0.5);

% --- Lyapunov Method ---
% Compute the matrix L representing the solution to the discrete Lyapunov equation
L_matrix = compute_lyapunov(a, b, g, r, T, Nh, Nk);

% Compute the functional value Phi using the Lyapunov solution matrix
Phi_Lyapunov = compute_lyapunov_functional(L_matrix, x0);

fprintf('Approximation using Lyapunov method: Phi_L = %f\n', Phi_Lyapunov);

% --- Monte Carlo Method ---
% Compute the functional value Phi using Monte Carlo simulation
Phi_MC = mc(a, b, g, r, T, Nh, Nk, N_mc, x0);

fprintf('Approximation using Monte Carlo method: Phi_MC = %f\n', Phi_MC);
