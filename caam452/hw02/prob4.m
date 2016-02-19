h=(1/2).^[2:5];
err0_a = zeros(length(h),1);
err1_b = zeros(length(h),1);
err0_a = zeros(length(h),1);
err1_b = zeros(length(h),1);

for j=1:length(h)

    N=(1-h(j))/h(j);
    A=(1/h(j))*(2*eye(N,N)-diag(ones(N-1,1),1)-diag(ones(N-1,1),-1));
    xInt=h(j)*[1 : N]';
    xFull=h(j)*[0 : (N+1)]';

    %when the exact solution is u(x)=x*(1-x), 
    %then f(x) is 2 and u(0)=0, u(1)=0
    b_a= h(j)*2* ones(N,1); 
    b_a(1)=b_a(1)+h(j)*0; b_a(N)=b_a(N)+h(j)*0;

    %when the exact solution is u(x)=x*(1-x)*e^(-x^2), 
    %then f(x) = (1-2x-2x^2+2x^3)*e^(-x^2) and u(0)=0, u(1)=0
    b_b=2*h(j)*(2*xInt.^4-2*xInt.^3-5*xInt.^2+3*xInt+1).*exp(-(xInt.^2));
    b_b(1)=b_b(1)+h(j)*0; b_b(N)=b_b(N)+h(j)*0;
    
    %approximate solutions
    V_a_approx = A\(b_a);
    U_a_approx = [0; V_a_approx ; 0];
    U_a_exact = xFull.*(1-xFull);
    V_b_approx = A\(b_b);
    U_b_approx = [0; V_b_approx ; 0];
    U_b_exact = xFull.*(1-xFull).*exp(-xFull.^2);
    
    %L2 error (err0) and Enery error (err1) of approximations
    err0_a(j) = ((h(j)/2)* ((U_a_approx - U_a_exact)'.^2)*([1 ; 2*ones(N,1) ;1]))^(1/2);
    aux_a = ((U_a_approx - U_a_exact)')*diag(ones(N,1),2)*(U_a_approx - U_a_exact);
    err1_a(j) = (1/(8*h(j))* (((((U_a_approx - U_a_exact)').^2)*([2 ; 3 ; 4*ones(N-2,1) ; 3 ; 2]))-4*aux_a))^(1/2);
    err0_b(j) = ((h(j)/2)* ((U_b_approx - U_b_exact)'.^2)*([1 ; 2*ones(N,1) ;1]))^(1/2);
    aux_b = ((U_b_approx - U_b_exact)')*diag(ones(N,1),2)*(U_b_approx - U_b_exact);
    err1_b(j) = (1/(8*h(j))* (((((U_b_approx - U_b_exact)').^2)*([2 ; 3 ; 4*ones(N-2,1) ; 3 ; 2]))-4*aux_b))^(1/2);


    
    %plots of approximate solutions when h=h(j)
    subplot(length(h),2,2*(j-1)+1);
        plot(xFull,U_a_approx); 
        xlabel('Nodes'); ylabel('Approx u(Nodes)'); axis([0 1 0 .3]);
        title({'approx soln to -u"=f when u(x)=x*(1-x)';sprintf('with h=%f',h(j))});
        hold on;
    subplot(length(h),2,2*(j-1)+2);
        plot(xFull,U_b_approx);     
        xlabel('Nodes'); ylabel('Approx u(Nodes)'); axis([0 1 0 .3]);
        title({'approx soln to -u"=f when u(x)=x*(1-x)*e^{-x^2}';sprintf('with h=%f',h(j))});
        hold on;
        
    %plot exact solutions as dotted curves    
    xAll = h(length(h))*[0 : (1-h(length(h)))/h(length(h))]';
    subplot(length(h),2,2*(j-1)+1);
        plot(xAll,xAll.*(1-xAll),':'); 
        axis([0 1 0 .3]);
        hold on;
    subplot(length(h),2,2*(j-1)+2);
        plot(xAll,xAll.*(1-xAll).*exp(-xAll.^2),':'); 
        axis([0 1 0 .3]); 
        hold on;
        
    clear V_a_approx U_a_approx U_a_exact aux_a
    clear V_b_approx U_b_approx U_b_exact aux_b
    clear A b_a b_b xInt N xFull


end

figure

%generation of L^2 error (err0) and Energy norm error (err1)
subplot(2,2,1)
    plot([h],[err0_a])    
    hold on;
    plot([h],[err1_a],':')
    xlabel('h val'); ylabel('error')
    title('L^2 (solid) and Energy (dotted) norm errors versus h')
    axes
subplot(2,2,2)
    plot([h],[err0_b])
	hold on;
    plot([h],[err1_b],':')
    xlabel('h val'); ylabel('error')
    title('L^2 (solid) and Energy (dotted) norm errors versus h')

%generation of loglog plots of h versus err0 and err1
subplot(2,2,3)
    loglog(h(1:length(h)),err0_a); 
    hold on;
    loglog(h(1:length(h)),err1_a,':')
    loglog(h(1:length(h)),h(1:length(h)).^2,'--')
    xlabel('log(h_k)'); ylabel('log(error)');
    title({'loglog of h versus L^2 (solid) and Energy (dotted) norm errors';
        'with loglog of h versus h^2 (dashed) for reference'});
subplot(2,2,4)
    loglog(h(1:length(h)),err0_b); 
    hold on;
    loglog(h(1:length(h)),err1_b,':')
    loglog(h(1:length(h)),h(1:length(h)).^2,'--')
    xlabel('log(h_k)'); ylabel('log(error)');
    title({'loglog of h versus L^2 (solid) and Energy (dotted) norm errors';
        'with loglog of h versus h^2 (dashed) for reference'});

% commented out in printout
% %generation of numerical rates of convergence for err0 and err1
% subplot(3,2,5)
%     rates0_a = log(err0_a./circshift(err0_a,[1 -1]))/log(2);
%     rates1_a = log(err1_a./circshift(err1_a,[1 -1]))/log(2);
%     plot(h(2:length(h)),rates0_a(1:(length(h)-1)));
%     hold on;
%     plot(h(2:length(h)),rates1_a(1:(length(h)-1)),':');
%     xlabel('h_k'); ylabel('rate(h_k) = log(h_{K-1}/h_{k})/log(2)');
%     title('L^2 (solid) and Energy (dotted) norm rate of convergence');
% subplot(3,2,6)
%     rates0_b = log(err0_b./circshift(err0_b,[1 -1]))/log(2);
%     rates1_b = log(err1_b./circshift(err1_b,[1 -1]))/log(2);
%     plot(h(2:length(h)),rates0_b(1:(length(h)-1)));
%     hold on;
%     plot(h(2:length(h)),rates1_b(1:(length(h)-1)),':');
%     xlabel('h_k'); ylabel('rate(h_k) = log(h_{K-1}/h_{k})/log(2)');
%     title('L^2 (solid) and Energy (dotted) norm rate of convergence');
%     

    
    clear all
