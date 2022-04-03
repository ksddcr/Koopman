%%B1=B(7,:);
B1=BJ;
n2=length(B1)-1;
p=2;
for j=1:5
X=B1(:,1:end-1);
Y=B1(:,2:end);
G = zeros(p+1,p+1);
A = zeros(p+1,p+1);
psiX = zeros(n2,p+1);  %% observable for X
psiY = zeros(n2,p+1);  %% observable for Y
for i = 1:n2
    psiX(i,:) = basis(X(i),p);
    psiY(i,:) = basis(Y(i),p);
    G = G + psiX(i,:)' * psiX(i,:);
    A = A + psiX(i,:)' * psiY(i,:);
end
G = G / n2;
A = A / n2;
K = pinv(G) * A;
B2=zeros(1,p+1);
B2(1,2)=1;
C=K*B2';
D1(1,j)=psiY(end,:)*C;
D1=round(D1);
B1=[B1,D1(1,j)];
n2=n2+1;
end