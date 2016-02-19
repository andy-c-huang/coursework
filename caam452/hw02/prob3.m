delta = 10.^(-[0:3]);
N=49;
h=1/(N+1);

for j=1:length(delta)
    L=-delta(j)/(h^2)*(-2*eye(N,N)+diag(ones(N-1,1),1)+diag(ones(N-1,1),-1));
    D=1/(2*h)*(diag(ones(N-1,1),1)-diag(ones(N-1,1),-1));
    A=L+D;
    b=zeros(N,1);
    b(N)=delta(j)/(h^2)-1/(2*h);
    V=A\b;
    U=[0; V ; 1];

    %generation of approx solution to -eps*u"+u'=0 and u(0)=0, u(1)=1
    subplot(length(delta),2,2*(j-1)+1)
    plot(h*[0 : N+1],U);
    xlabel('Nodes'); ylabel('Approx u(Nodes)');
    hold on;
    title({'approx soln to -\epsilon *u"+u''=0 and u(0)=0, u(1)=1';
            sprintf('when \\epsilon=%f',delta(j))})
    plot(h*[0 : N+1], (exp((h*[0 : N+1])*(1/delta(j)))-1)/(exp(1/delta(j))-1),':')
    %generation of exact solution
    subplot(length(delta),2,2*(j-1)+2)
    plot(h*[0 : N+1], (exp((h*[0 : N+1])*(1/delta(j)))-1)/(exp(1/delta(j))-1))
    xlabel('Nodes'); ylabel('Exact u(Nodes)');
    title({'exact soln to -\epsilon *u"+u''=0 and u(0)=0, u(1)=1';
            sprintf('when \\epsilon=%f',delta(j))})
    
    clear L D A b U
    
end
    