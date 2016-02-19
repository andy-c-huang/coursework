function [S]=hw04prob3b(E)
S = zeros(3,3); % stiffness matrix

A1=E(1,:);
A2=E(2,:);
A3=E(3,:);
B=[(A2(1)-A1(1)), (A3(1)-A1(1)); (A2(2)-A1(2)), (A3(2)-A1(2))];
delta = det(B);
a = B(1,1);
b = B(1,2);
c = B(2,1);
d = B(2,2);

DPhi = [-1, -1; 1, 0; 0, 1];
% DPhi is the total derivative of Phi = (Phi^ _1, Phi^ _2, Phi^ _3)
% where Phi^ _k = 1-s-t,s,t for k=1,2,3
% e.g., DPhi(2,1) = del(Phi^ _2)/del(t)

for i = 1:3
    for j = 1:3
        S(i,j) = (1/(2*delta))*((d^2+b^2)*DPhi(i,1)*DPhi(j,1)+(-d*c-a*b)*DPhi(i,1)*DPhi(j,2)+(-d*c-a*b)*DPhi(i,2)*DPhi(j,2) + (c^2+a^2)*DPhi(i,2)*DPhi(j,2) );
        % using the expansion from (a), we have the above line equal to
        % S(i,j) = (delta/2)*((inv(B)')*((DPhi(i,:))'))'*((inv(B)')*((DPhi(j,:))'));
    end
end