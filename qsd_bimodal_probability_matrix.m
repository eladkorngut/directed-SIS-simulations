function [tau,prob_matrix] = qsd_bimodal_probability_matrix(N,R0,gamma,epsilon_lam,epsilon_mu,k)
%  Caluclate the master eq for 2d sis model, in case of large N
%   Returns tau which is the extinction time and prob_matrix which is the
%   quasi-stationry distrubition

mu1=1-epsilon_mu;
mu2=1+epsilon_mu;
lam1=1-epsilon_lam;
lam2=1+epsilon_lam;
beta=R0/(k);
beta=beta/(1+epsilon_lam*epsilon_mu);
L=((N/2+1))^2;
i=[];
j=[];
v=[];
for col=2:L-1
    c1=col;
    c2=mod(col+(N/2)+1,L+1);
    c3=mod(col-(N/2)-1,L);
    I1=double(idivide(col-1,int16((N/2+1))));
    I2=mod(col-1,N/2+1);
    i(end+1)=c1;
    j(end+1)=col;
    v(end+1)= -(beta*(lam1*I1+lam2*I2)*(mu1*(N/2-I1)+mu2*(N/2-I2))+gamma*I1+gamma*I2);
    if (I2+1<=N/2)
        i(end+1)=c1+1;
        j(end+1)=col;
        v(end+1)=gamma*(I2+1);
    end
    if(I1+1<=N/2)
        i(end+1)=c2;
        j(end+1)=col;
        v(end+1)=gamma*(I1+1);
    end
    if(I1-1>=0 && I2>=0)
        i(end+1)=c3;
        j(end+1)=col;
        v(end+1)=beta*(lam1*(I1-1)+lam2*I2)*mu1*(N/2-(I1-1));
    end
    if(I2-1>=0 && I1>=0)
        i(end+1)=c1-1;
        j(end+1)=col;
        v(end+1)=beta*(lam1*I1+lam2*(I2-1))*mu2*(N/2-(I2-1));
    end
end
i(end+1)=L;
j(end+1)=L;
v(end+1)=-(gamma*(N/2)+gamma*(N/2));
i(end+1)=L-1;
j(end+1)=L;
v(end+1)=beta*(lam1*N/2+lam2*(N/2-1))*mu2*(N/2-(N/2-1));
i(end+1)=L-N/2-1;
j(end+1)=L;
v(end+1)=beta*(lam1*(N/2-1)+lam2*N/2)*mu1*(N/2-(N/2-1));
Q = sparse(i,j,v,L,L);
Q(1,:)=[]; Q(:,1)=[];
Q=transpose(Q);
[prob_vec,D] = eigs(Q,1,'smallestabs');
tau = -1/D;
prob_vec =prob_vec/sum(prob_vec);
prob_vec=cat(1,[0],prob_vec);
for col=1:L
    I1=double(idivide(col-1,int16((N/2+1))))+1;
    I2=mod(col-1,N/2+1)+1;
    prob_matrix(I1,I2)=prob_vec(col);
end

end