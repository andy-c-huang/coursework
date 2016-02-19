function [ maxnormerror ] = prob4( M )
c=[0,1,1];
for m=1:length(M)
   clear h U A a nodes b
   h=1/M(m);
   U=zeros(M(m)-1,3);
   A=zeros(3,M(m)-1,M(m)-1);
   nodes = (1:(M(m)-1))*h;
   b(1,:)=  (0*nodes); b(1,1)=b(1,1)-1*(1/h^2); b(1,M(m)-1)=b(1,M(m)-1)-0*(1/h^2);
   b(2,:) = (1*nodes); b(2,1)=b(2,1)-0*(1/h^2); b(2,M(m)-1)=b(1,M(m)-1)-1*(1/h^2);
   b(3,:) = (2*(nodes-2).*exp(-nodes)); b(3,1)=b(3,1)-(0*exp(-0))*(1/h^2); b(3,M(m)-1)=b(3,M(m)-1)-(1*exp(-1))*(1/h^2);
      %these are the appropriate functions f(x) evaluated at the M-1 points: 1*h, 2*h, \ldots, (M-1)*h
      %here, we check the FDA technique against known exact solutions 
      %u(x) = 1-x and c=0 (so, u''+0*u = f(x)=0)
      %u(x)=x and c=1 (so, u''+1*u = f(x)=x)
      %u(x)=x*e^(-x) and c=1 (so, u''+1*u = f(x)=f(x)=2*(x-1)*e^(-x)) 
   for k=1:length(c)
      h=1/M(m);
      A(k,:,:) = c(k)*eye(M(m)-1,M(m)-1) + 1/h^2 *(diag(-2*ones(M(m)-1,1),0)+diag(ones(M(m)-2,1),1)+diag(ones(M(m)-2,1),-1));
         %there are M-1 entries because, among u(0*h), u(1*h), \ldots, u((M-1)*h), u(M*h=1), we don't record u(0) and u(1)
      a(:,:)=A(k,:,:);
      U(:,k)=a\(b(k,:)');
   end
maxnormerror(1,m) = max(abs(U(:,1)-(1-nodes)'));
maxnormerror(2,m) = max(abs(U(:,2)-(nodes)'));
maxnormerror(3,m) = max(abs(U(:,3)-(nodes.*exp(-nodes))'));
end

rates = log(maxnormerror(3,:)./circshift(maxnormerror(3,:),[1 -1]))/log(2);
   %make sure to delete the last entry of this rates vector, as it contains bogus information
steps = 1./M;
plot(steps(2:length(M)),rates(1:(length(M)-1)))
   %make sure to omit the first rate of 1/10

end

