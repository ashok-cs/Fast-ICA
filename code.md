```matlab
% MATLAB CODE FOR FAST ICA ALGORITHM
%epsilon -> convergence factor
epsilon=0.01;   
Niter= 100;

%load input signals data X
load X; 
[M N] =size(X);
B = zeros(M);
C = cov (X');
[V, D] = eig(C);
% V    matrix of eigenvectors of C
% D    Diagonal matrix of eigenvalues of C
Q = inv(sqrt(D))*V'; % Whitening Matrix
Qinv = V*sqrt (D);
At=Q*A;
Xt=Q*X;
%------FAST ICA ITERATION UNIT---------------------
for m=1:M;
    %initial Vector
    w = randn(M,1); w=w-B*B'*w;  w=w/norm(w);
    w_old = zeros(size(w));
    for iter=1:Niter, 
      Xg=[]; g_der_vec=[]; Gvec=[];
      w=w-B*B'*w; w=w/norm(w);
      % Test for convergence
      if norm(w - w_old)<epsilon | norm(w + w_old)<epsilon,
      B(:,m)=w; Wt(m,:)=w'*Q; Wtinv(:,m)=Qinv*w;
      break
      end
           for j=1:N;
           x=Xt(:,j);  
           u=x'*w; G= -exp(-u^2/2); g=u*exp(-u^2/2);
           g_der=(1-u^2)*exp(-u^2/2); Xg=[Xg x*g];
           g_der_vec=[g_der_vec   g_der]; Gvec= [Gvec G];
           end;
      w_old=w;
      w=mean(Xg')'-(mean(g_der_vec)*w); GG(iter)=mean(Gvec);Gsave=Gvec;
      w=w/norm(w);
    end
 end

disp('Wtinv*Wt');Wtinv*Wt
Shat=Wt*X;
Shat1=Shat(1,:); Shat2=Shat(2,:);if M==3, Shat3=Shat(3,:); end;

% shat -> estimated source signal
% -------------------------------------------------------------
```
