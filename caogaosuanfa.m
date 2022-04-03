for i=1:10
    B2(i,:)=D3(1,i:end-10+i);
end
X = B2(:,1:end-1);
Y= B2(:,2:end);%%B为输入原始矩阵，X为系统前n-1个时刻状态，Y为系统后n-1个时刻的状态
for j=1:5
[U, S, V] = svd( X, 'econ' );%%对X进行简易的奇异值分解
E= pinv(S);%%计算X的奇异值分解中S矩阵的广义逆矩阵
Aw = U'*Y*V*E;%%定义的A弯
[w, u] = eigs(Aw, size(Aw,1));%%计算A弯的特征值与特征向量
Modes=U*w;%%得出DMD的模态
A=Modes*u*(pinv(w))*U';%%相应的DMD近似的Koopman算子可以得出
C=A*Y;
C=real(C);
D=round(C);
%%F(:,j)=sum(D(:,end))/2;
F(:,j)=D(end,end);%%一步预测，取实数部分，取整。
X=[X,Y(:,end)];
Y=[Y,D(:,end)];
end


for i=1:2
    for j=1:9
        B1(i,j)=BJ(1,end-(i-1)-(9-j)*i);
    end
end%%第二种升维方法



D1=D(:,1);
D2=D(:,2);
D3=D2-D1;
D3=abs(D3);
D4=D3./D2;
D5=sum(D4);
D6=D5/20;%%MAPE指标

D1=D(:,1);
D2=D(:,2);
D3=D2-D1;
D3=abs(D3);
D4=D3.^2;
D5=sum(D4);
D6=D5/20;
D7=D6^(1/2);%%RMSE指标




%%B1=B(7,:);
for k=1:34
B1=B(k,:);
n2=length(BJ)-1;
p=2;
for j=1:5
X=BJ(:,1:end-1);
Y=BJ(:,2:end);
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
D1(k,j)=psiY(end,:)*C;
D1=round(D1);
B1=[B1,D1(k,j)];
n2=n2+1;
end
B1=zeros
end

D2=D(1,1:end-1);
D3=D(1,2:end);
D4=D3-D2;
D5=D4(1,end-59:end);

x = linspace(datenum('2020/03/01'),datenum('2021/03/31'),396);   %生成376个时间
y = D;    % 生成时间序列y
plot(x,y);          % 画图
set(gca,'XTick',...
    [737851  737879  737907  737935  737963  737991  738019  738047  738075  738103  738131  738159  738187  738215  738243],...
    'XTickLabel',...
    {'03-01','03-29','04-26','05-24','06-21','07-19','08-16','09-13','10-11','11-08','12-06','01-03','01-31','02-28','03-28'});
plot(D2);
set(gca,'XTick',...
    [1 5 10 15 20 ],...
    'XTickLabel',...
    {'06-11','06-15','06-20','06-25','06-30'});


set(gca,'XTick',...
    [1 57 113 169 225 281 337 393],...
    'XTickLabel',...
    {'04-01','05-27','07-22','09-16','11-11','01-06','03-03','04-28'});





plot(D);
set(gca,'XTick',...
    [1 85 169 253 337 421],...
    'XTickLabel',...
    {'04-01','06-24','09-16','12-09','03-03','05-26'});

D1=2161591
for i=1:20
    D1=D1+D(1,i);
    D2(1,i)=D1;
end


plot(D);
set(gca,'XTick',...
    [1 57 113 169 225 281 337 393],...
    'XTickLabel',...
    {'05-01','06-26','08-21','10-16','12-11','02-05','04-03','05-29'});



for j=1:407
     B=D(:,j:j+29);
    for i=1:10
     Q(i,:)=B(1,i:i+end-10);
    end
InputSnapshots = Q(:,1:end-1);
OutputSnapshots = Q(:,2:end);
[Ux, Sx, Vx] = svd( InputSnapshots, 'econ' );
SxInv = pinv(Sx);
Atilde = Ux' * OutputSnapshots * Vx * SxInv;
[w, lambdas] = eigs(Atilde, size(Atilde,1));
Modes = Ux * w;
x=diag(lambdas); %%K=Modes*lambdas*(pinv(w))*Ux';
INM=(pinv(w))*Ux';
for k=1:10
    y1(:,k)=norm(Modes(:,k));
end
all_y=sum(y1);
for l=1:10
    y2(:,l)=y1(:,l)/all_y;
end
y3(:,j)=y2';
y4(:,j)=x;
end
y5=real(y4);
for p=1:407
    y7(1,p)=max(y5(:,p));
end

for i=1:14
    y1(i,1)=1/14;
end

for i=1:406
    f1=x2(:,i);%%这里是第一个地区特征值的矩阵
    f2=x4(:,i);%%这里是第二个地区特征值的矩阵
    w1=y1;%%这里是第一个地区的特征向量算得的权重
    w2=y4;%%这里是第二个地区的特征向量算得的权重
    [x fval] = emd_KEM(f1, f2, w1, w2, @gdf);
    Dis(:,i)=fval;
    Ju(:,i)=x;
end