clc; clear; format Long G

G=pi/200;
D=[104  102 666.9285
   106  102 1064.68565
   105  104 627.647
   106  104 750.3989
   106  105 326.3115
   102  105 1136.005];
D=sortrows(sortrows(D,2),1,"descend");

X=[105 31.8140
   104 23.46216
   102 -46.0387];
B106=873.7409;

x=[B106+X(:,2)]';

l=[106 105   31.8140
   106 104   23.46216
   106 102   -46.0387
   102 105   77.8395
   104 102   -69.54468
   105 104   -8.31098];
l1=sortrows(sortrows(l,2),1,"descend");

ST1="B"+num2str(l1(:,1));
CL1="B"+num2str(l1(:,2));
AR1="B"+num2str(sortrows(X(:,1),1,"descend"));
P=eye(size(l,1)).*(1./(((D(:,3))/1000).^2));

l=l1(:,3);

AR=str2sym(AR1);
ST=str2sym(ST1);
CL=str2sym(CL1);
for n=1:size(ST,1)
    for m=1:size(x,2)
        A(n,m)=diff(ST(n)-CL(n),AR(m));
    end
end
A=eval(A);

l=l*(-1);

xA=(A'*P*A)^(-1)*A'*P*l;
v=A*xA-l;
lc=l+v;

Vys_vysky=B106+xA;

s0=sqrt((v'*P*v)/3);
Qx=(s0^2)*inv(A'*P*A);       %kovarianční matice neznámých
Qh=(s0^2)*A*inv(A'*P*A)*A';  %kovarianční matice vyrovnaných převýšení

SM_nez_X=sqrt(diag(Qx))';     
SM_vyr_mer=sqrt(diag(Qh))';

alf1=0.05;                  alf2=0.95;
kv1=chi2inv(alf1,3);        kv2=chi2inv(alf2,3);
mez1=sqrt(kv1/3);           mez2=sqrt(kv2/3);
