
% solve the Poisson problem u_{xx} + u_{yy} = f(x,y)
% on [0,1] x [0,1].  
% 
% The 5-point Laplacian is used at interior grid points.
% 

%has been modified from script to take any input value of N
%h = 1/10,1/20,,1/40,1/80
function[err]= poisson(N)
%N = 9;
h = 1/(N+1);
x = linspace(0,1,N+2);   % grid points x including boundaries
y = linspace(0,1,N+2);   % grid points y including boundaries


[X,Y] = meshgrid(x,y);      % 2d arrays of x,y values
X = X';                     % transpose so that X(i,j),Y(i,j) are
Y = Y';                     % coordinates of (i,j) point

Iint = 2:N+1;              % indices of interior points in x
Jint = 2:N+1;              % indices of interior points in y
Xint = X(Iint,Jint);       % interior points
Yint = Y(Iint,Jint);

f = @(x,y) 1.25*exp(x+y/2);         % f(x,y) function

rhs = f(Xint,Yint);        % evaluate f at interior points for right hand side
                           % rhs is modified below for boundary conditions.

utrue = exp(X+Y/2);        % true solution for test problem

% set boundary conditions around edges of usoln array:

usoln = utrue;              % use true solution for this test problem
                            % This sets full array, but only boundary values
                            % are used below.  For a problem where utrue
                            % is not known, would have to set each edge of
                            % usoln to the desired Dirichlet boundary values.


% adjust the rhs to include boundary terms:
rhs(:,1) = rhs(:,1) - usoln(Iint,1)/h^2;
rhs(:,N) = rhs(:,N) - usoln(Iint,N+2)/h^2;
rhs(1,:) = rhs(1,:) - usoln(1,Jint)/h^2;
rhs(N,:) = rhs(N,:) - usoln(N+2,Jint)/h^2;


% convert the 2d grid function rhs into a column vector for rhs of system:
F = reshape(rhs,N*N,1);

% form matrix A:
I = speye(N);
e = ones(N,1);
T = spdiags([e -4*e e],[-1 0 1],N,N);
S = spdiags([e e],[-1 1],N,N);
A = (kron(I,T) + kron(S,I)) / h^2;

% Solve the linear system:
uvec = A\F;  

% reshape vector solution uvec as a grid function and 
% insert this interior solution into usoln for plotting purposes:
% (recall boundary conditions in usoln are already set) 

usoln(Iint,Jint) = reshape(uvec,N,N);

% assuming true solution is known and stored in utrue:
err = max(max(abs(usoln-utrue)));   
fprintf('Error relative to true solution of PDE = %10.3e \n',err)

% plot results:

clf
hold on

% plot grid:
plot(X,Y,'k');  plot(X',Y','k')

% plot solution:
contour(X,Y,usoln,30)

axis([0 1 0 1])
daspect([1 1 1])
title('Contour plot of computed solution')
hold off

figure
title('Computed solution')
surf(X,Y,usoln)

