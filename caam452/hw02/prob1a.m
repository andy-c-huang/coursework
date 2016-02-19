N=[9; 19; 39; 79]
err=zeros(length(N),1);

for i=1:length(N)
    h = 1/(N(i)+1);
    x = linspace(0,1,N(i)+2);   % grid points x including boundaries
    y = linspace(0,1,N(i)+2);   % grid points y including boundaries


    [X,Y] = meshgrid(x,y);      % 2d arrays of x,y values
    X = X';                     % transpose so that X(i,j),Y(i,j) are
    Y = Y';                     % coordinates of (i,j) point

    Iint = 2:N(i)+1;              % indices of interior points in x
    Jint = 2:N(i)+1;              % indices of interior points in y
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
    rhs(:,N(i)) = rhs(:,N(i)) - usoln(Iint,N(i)+2)/h^2;
    rhs(1,:) = rhs(1,:) - usoln(1,Jint)/h^2;
    rhs(N(i),:) = rhs(N(i),:) - usoln(N(i)+2,Jint)/h^2;


    % convert the 2d grid function rhs into a column vector for rhs of system:
    F = reshape(rhs,N(i)*N(i),1);

    % form matrix A:
    I = speye(N(i));
    e = ones(N(i),1);
    T = spdiags([e -4*e e],[-1 0 1],N(i),N(i));
    S = spdiags([e e],[-1 1],N(i),N(i));
    A = (kron(I,T) + kron(S,I)) / h^2;

    % Solve the linear system:
    uvec = A\F;  

    % reshape vector solution uvec as a grid function and 
    % insert this interior solution into usoln for plotting purposes:
    % (recall boundary conditions in usoln are already set) 

    usoln(Iint,Jint) = reshape(uvec,N(i),N(i));

    % assuming true solution is known and stored in utrue:
    err(i) = max(max(abs(usoln-utrue)));   
    fprintf('Error relative to true solution of PDE = %10.3e \n',err(i))

    % plot grid:
    subplot(length(N),2,2*(i-1)+1);
    plot(X,Y,'k');  plot(X',Y','k')

    % plot solution:
    contour(X,Y,usoln,30)

    axis([0 1 0 1])
    daspect([1 1 1])
    title(sprintf('Contour plot of computed solution when h=%f',h))
    hold on

    subplot(length(N),2,2*(i-1)+2);
    title(sprintf('Computed solution when h=%f',h))
    surf(X,Y,usoln)
    hold on
    
    clear h x y X Y Iint Jint Xint Yint rhs utrue usoln F I e T S A uvec
end

figure
subplot(2,1,1);
    plot(1./(N+1),err)
    xlabel('h val'); ylabel('error')
    title('abs max norm error = |u_{exact}-u_{approx}|')
    strValues = strtrim(cellstr(num2str([(1./(N+1)) err],'(%d,%d)')));
    text((1./(N+1)), err,strValues,'VerticalAlignment','bottom');
subplot(2,1,2);
    loglog(1./(N+1),err);
    hold on;
    loglog(1./(N+1),(1./(N+1)).^2,':');
    xlabel('log(hval)'); ylabel('log(abs max norm error)')
    title({'(loglog of hval versus max norm error (solid)'; 'loglog of hval versus hval^2 (dotted) for reference)'})
    
    
%     rates= log(err./circshift(err,[1 -1]))/log(2);
%     plot(N(2:length(N)),rates(1:(length(N)-1)));