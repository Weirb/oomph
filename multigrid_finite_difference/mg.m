function [] = gmg_final()


    % Plot the solution/iteration data
    plots = true;


    % Domain 
    global a b; 
    a = 0.1; b = a+1;


    % Global functions for mg. Conditionally set these based on problem.
    global A create_mat smooth;
    
    % Grid resolution
    k = 10;
    N = 2^k + 1;
    
    % Grid width
    h = (b-a)/(N+1);
    
    % Interval data
    x = linspace(a,b,N)';
    
    % Variable coefficient or Poisson?
    poisson = 0;


    % If we want to solve Poisson problem
    if poisson
        % Boundary conditions
        alpha = 0; beta = 1;


        % Exact solution
        u_ex = (alpha-beta)/(a-b)*x + 0.5*(alpha+beta) + 0.5*(a+b)/(a-b)*(beta-alpha);
        f = 0*x;
        f(1) = f(1) - alpha/h^2;
        f(N) = f(N) - beta/h^2;


        A = @(s) A_poisson(s);
        create_mat = @(s) create_mat_poisson(s);
        smooth = @(u,v,w) smooth_poisson(u,v,w);
     
    else % Otherwise, we are solving variable coefficient problem


        % Exact solution
        u_ex = log(x);


        % Boundary conditions
        alpha = log(a); beta = log(b);
        f = 0*x;
        f(1) = f(1) - alpha/h^2;
        f(N) = f(N) - beta*(1/h^2+1/(a+N*h)/h);
        
        A = @(s) A_varcoef(s);
        create_mat = @(s) create_mat_varcoef(s);
        smooth = @(u,v,w) smooth_varcoef(u,v,w);
    end
    
    % Tolerances
    res_tol = 1e-10;
    ratio_tol = 1e-2;
    
    % Arrays to hold convergence data
    residuals = [];
    errors = [];
    
    % 1 for V-cycle, 2 for W-cycle
    v_or_w = 1;
    
    % How many iterations?
    max_cycles = 100;
    
    % Begin mg iteration
    u = zeros(N,1);
    for cycles = 1:max_cycles
        
        % Perform a single mg cycle
        u = mg_cycle(u, f, k, v_or_w);


        % Store the current residual and error
        residuals = [residuals, norm(A(u)-f)];
        errors = [errors, norm(u-u_ex)];
        
        % Check if we have reached convergence tolerance
        if residuals(end) < res_tol
            fprintf('\n  **CONVERGED**\n\n');
            break;
        end
        
        % Check if we are no longer reducing the residual
        if cycles > 1
            ratio = residuals(cycles)/residuals(cycles-1);
            if norm(ratio-1,2) < ratio_tol
                fprintf('\n  **RESIDUAL NOT REDUCING**\n\n');
                break;
            end
        end
        
        % Check if something has gone terribly wrong
        if any(isnan(residuals))
            fprintf('\n  **FAILED**\n\n');
            break;
        end
    end
    
    % Display convergence information
    if v_or_w == 1; c = 'V'; else c = 'W'; end;
    fprintf('\t%c-cycle count: %d\n',c,cycles);
    fprintf('\tFinal residual: %.3e\n',residuals(end));
    fprintf('\tFinal error   : %.3e\n\n',errors(end));
    
    if plots
        % Plot the numerical and exact solutions
        figure();
        subplot(211);
        hold on
        plot(x, u, 'k.')
        plot(x, u_ex,'r','LineWidth',2);
        xlabel('x');
        ylabel('u(x)');
        title('SOLUTION');
        legend('Calculated', 'Exact','Location','SouthEast');


        % Plot the residual and error against iteration count
        subplot(212);
        set(gca,'YScale','log')
        hold on;
        plot(residuals,'ko-');
        plot(errors,'ro-');
        xlabel('Iteration number');
        ylabel('Residual');
        title('ITERATION COUNT');
        legend('Residual', 'Error');
        grid on;
    end
end


function u = mg_two_grid(u,f)
    global A create_mat smooth;


    % Pre-smooth
    u_fine = smooth(u, f, jacobi_pre);
    
    % Restrict the residual onto the coarser level
    residual = restrict(A(u_fine) - f);
    
    % Solve the problem exactly on the coarse level
    z = create_mat(2^(n-1)+1);
    u_coarse = z\residual;
    
    % Apply the coarse grid correction
    u_fine = u_fine - prolong(u_coarse);
    
    % Post-smooth
    u = smooth(u_fine, f, jacobi_post);
end


function u = mg_cycle(u, f, n, w)
    global A create_mat smooth;


    % w = 1 gives the v-cycle
    % w = 2 gives the w-cycle
    
    % Number of smoothing iterations
    jacobi_pre = 100;
    jacobi_post = 100;
    
    % Base case, solve the problem exactly
    if n == 3
        z = create_mat(2^n+1);
        u = z\f;
    else
        % Pre-smooth
        u = smooth(u, f, jacobi_pre);
        
        % Compute the residual and restrict to coarser grid
        d = restrict(A(u) - f);
        v = 0*d;
        for a = 1:w
            % Compute w recursive steps using the residual
            v = mg_cycle(v, d, n-1, w);
        end
        
        % Add the correction back to the finer grid solution
        u = u - prolong(v);
        
        % Post-smooth
        u = smooth(u, f, jacobi_post);
    end
end


function u = mg_f_cycle(f, n_max, n_min)
    % Perform a multigrid F-cycle, reaching the finest level n_max with 
    % 2^n_max+1 elements and starting on the coarsest level n_min


    % Solve exactly on the coarsest level first
    z = create_mat(2^n_min+1);
    u = z\restrict_n(f,n_max-n_min);
    
    % Perform n V-cycles, moving up a level each iteration
    for i = n_min:n_max
        % Interpolate to finer mesh first
        u = prolong(u);
        % Perform several V-cycles
        for k = 1:20
            u = mg_cycle(u, restrict_n(f, n_max-i), i, 1);
        end
    end
end


function u = smooth_varcoef(u, f, max_iter)
% Smoothing operation is done by the weighted Jacobi iteration


    % Variables a and b are the domain end points
    global a b
    % Weight parameter
    w = 2/3;
    % Number of elements on the current level
    m = numel(u);
    % Grid size on the current level
    h = (b-a)/(m+1);
    % Array containing all of the interior points
    I = 2:m-1;
    % Variable for tracking old solution
    u_old = u;


    % Iterate until reached maximum number of iterations
    for k=1:max_iter
        % Perform weighted Jacobi iteration
        % Treat end points separately
        % Vectorised operation to improve performance
        u(1) = (1-w)*u(1) - w/(-2/h^2 - 1/(a+h)/h ).*(u(2)*(1/h^2+1/(a+h)/h) - f(1));
        u(I) = (1-w)*u(I) - w./(-2/h^2 - 1./(a+I'*h)/h ).*(u(I-1)/h^2 + u(I+1).*(1/h^2+1./(a+I'*h)/h) - f(I));
        u(m) = (1-w)*u(m) - w/(-2/h^2 - 1/(a+m*h)/h ).*(u(m-1)/h^2 - f(m));
        
        % If we are no longer reducing the solution, stop iteration
        if norm(u_old - u) < 1e-10
            break;
        end
        % Keep the current solution for next time
        u_old = u;
    end
end


function u = smooth_poisson(u, f, max_iter)
% Smoothing operation is done by the weighted Jacobi iteration


    % Variables a and b are the domain end points
    global a b
    % Weight parameter
    w = 2/3;
    % Number of elements on the current level
    m = numel(u);
    % Grid size on the current level
    h = (b-a)/(m+1);
    % Array containing all of the interior points
    I = 2:m-1;
    % Variable for tracking old solution
    u_old = u;
    
    % Iterate until reached maximum number of iterations
    for k=1:max_iter
        % Perform weighted Jacobi iteration
        % Treat end points separately
        % Vectorised operation to improve performance
        u(1) = (1-w)*u(1) + w/2*(u(2) - h^2*f(1));
        u(I) = (1-w)*u(I) + w/2*(u(I-1) + u(I+1) - h^2*f(I));
        u(m) = (1-w)*u(m) + w/2*(u(m-1) - h^2*f(m));
        
        % If we are no longer reducing the solution, stop iteration
        if norm(u_old - u) < 1e-10
            break;
        end
        % Keep the current solution for next time
        u_old = u;
    end
end


function v = restrict(v)
    % Restriction interpolates from the fine grid down to the coarse grid
    v = v(1:2:end);
end


function v = restrict_n(v, n)
    % Iteratively apply restriction to a vector to reach the desired level
    % Used in the mg F-cycle
    for i=1:n
        v = restrict(v);
    end
end


function x = prolong(v)
    % Prolongation interpolates from the coarse grid up to the fine grid
    n = numel(v);
    x = interp1(linspace(0,1,n)', v, linspace(0,1,n*2-1)', 'linear');
end


function x = A_varcoef(u)
% Compute the matrix action for the variable coefficient discretisation
    
    % Domain end points, a and b
    global a b
    % Number of elements on the current level
    N = numel(u);
    % Grid size on the current level
    h = (b-a)/(N+1);
    % Solution variable
    x = u*0;
    
    % Apply the discretisation stencil to the vector u
    % Treat end points separate from interior points
    x(1) = [-2/h^2 - 1/(a+h)/h, 1/h^2 + 1/(a+h)/h]*u(1:2);
    for i = 2:N-1
         x(i) = [1/h^2, -2/h^2-1/(a+i*h)/h, 1/h^2+1/(a+i*h)/h]*u(i-1:i+1);
    end
    x(N) = [1/h^2, -2/h^2 - 1/(a+N*h)/h]*u(N-1:N);
end


function A = create_mat_varcoef(N)
% Create matrix of size NxN for solving the variable coefficient problem exactly


    % Domain end points, a and b
    global a b
    % Grid size on level N
    h = (b-a)/(N+1);
    % Matrix variable
    A = zeros(N);
    
    % Store the matrix entries into A
    A(1,1:2) = [-2/h^2 - 1/(a+h)/h, 1/h^2 + 1/(a+h)/h];
    for j=2:N-1
        A(j,j-1:j+1) = [1/h^2 , -2/h^2 - 1/(a+j*h)/h, 1/h^2 + 1/(a+j*h)/h];
    end
    A(N,N-1:N) = [1/h^2, -2/h^2 - 1/(a+N*h)/h];
end


function x = A_poisson(u)
% Compute the matrix action for the Poisson discretisation


    global a b;
    N = numel(u);
    h = (b-a)/(N+1);
    x = 0*u;


    x(1) = 1/h^2 * (-2*u(1) + u(2));
    for i = 2:N-1
        x(i) = 1/h^2 * (u(i-1) -2*u(i) + u(i+1));
    end
    x(end) = 1/h^2 * (u(end-1) - 2*u(end));


end


function A = create_mat_poisson(n)
% Create matrix of size NxN for solving the Poisson problem exactly


    global a b;
    h = (b-a)/(n+1);
    A = 1/h^2 *(sparse(diag(-2*ones(n,1))) + sparse(diag(ones(n-1,1),1)) + sparse(diag(ones(n-1,1),-1)));
end