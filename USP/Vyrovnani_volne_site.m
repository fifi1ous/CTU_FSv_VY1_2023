clc; clear; format Long G

global G2R;
G2R=pi/200;
%% Inicializace hodnot
%   B           Y               X
SS=[102     845560.366     998311.676     
    104     845326.159     997687.739     
    105     845703.651     997183.642     
    106     845993.693     997336.573];

%       ST  CL  SM
smery=[102	106	0
       102	104	49.4886
       104	106	0
       104	105	28.2515
       104	102	292.0341
       105	104	0
       105	106	110.0354
       106	105	0
       106	104	61.7133
       106	102	104.2623];

delky=[102	104	666.459
       102	106	1067.067
       104	105	629.776
       104	106	754.320
       105	106 327.934];
zach=[ 105  102];
%% Nastavovací výpočty
SS=sortrows(SS,1);
smery=sortrows(sortrows(smery,3),1);
delky=sortrows(sortrows(delky,2),1);

smery(:,3)=smery(:,3)+5;
for n=1:size(smery,1)    
        st=find(smery(n,1)==SS(:,1));
        cl=find(smery(n,2)==SS(:,1));
        OP(n)=atan2(SS(cl,2)-SS(st,2),SS(cl,3)-SS(st,3))-smery(n,3)*G2R;
        if OP(n)<0
            OP(n)=OP(n)+2*pi;
        end
end
smery2=[smery,OP'];

OP=[];
for n=1:size(SS,1)
    st=find(SS(n,1)==smery(:,1));
    OP=[OP;mean(smery2(st,4))];
end

SS=[SS,OP];
[X] = Prevod_SS_M_na_V(SS);
%% výpočet vah
sm_u=0.6/1000*G2R;
sm_d=(1.5+2*(delky(:,3)/1000))/1000;
P1=eye(size(smery,1))*(1/(sm_u^2)); PN=zeros(size(smery,1),size(delky,1));
P2=eye(size(delky,1)).*(1./(sm_d.^2)); 
P=[P1,PN;PN',P2];
%% Vyrovnání
for n=1:10
    [lc,lp,lm] = vypocet_l(smery,delky,SS);
    [BT,b,b1] = Matice_B_b(zach,SS,X);
    [A] = MAT_A(X,SS,smery,delky);
    if n>2
        if v-lc<10e-9
            break
        end
    end
    ATPA=[A'*P*A,BT';BT,b];
    ATPL=[A'*P*lc;b1];
    dxk=ATPA^(-1)*ATPL;
    dx=dxk(1:end-3);
    v=A*dx+lc;
    X=X+dx';
    [SS] = Prevod_SS_V_na_M(X,SS);
end
SS(:,4)=SS(:,4)/G2R+5;
%% Přesnosti
s0=sqrt((v'*P*v)/(size(A,1)+size(BT,1)-size(A,2)));
Mx=diag(ATPA^(-1));
Mx=sqrt(Mx(1:end-3,:));
[Mx] = Prevod_SS_V_na_M(Mx'*1000,SS); Mx(:,4)=Mx(:,4)/G2R;
alf1=0.05;                  alf2=0.95;
kv1=chi2inv(alf1,3);        kv2=chi2inv(alf2,3);
mez1=sqrt(kv1/3);           mez2=sqrt(kv2/3);
%% Výpis výsledků
if s0>=mez1 && s0<=mez2
    fprintf("Gratuluji byl jsi vyrovnán jsi v zenu  (^‿^)\n");load handel.mat;nBits = 16;sound(y,Fs,nBits);
    fprintf("\nAposteriorní sm. odch. jednotková: %4.1f\n",s0);
    fprintf("%4.1f <= %4.1f <= %4.1f\n",[mez1,s0,mez2]);
    fprintf("\nSeznam vyrovnaných souřadnic:\n CB.    Y[m]         X[m]       OP[g]\n")
    fprintf("%d %12.3f %12.3f %8.4f\n",SS');
    fprintf("Seznam směrodatných odchylek vyrovnaných souřadnic:\n CB.  Y[mm]     X[mm]    OP[mgon]\n")
    fprintf("%d  %5.2f     %5.2f     %5.2f\n",Mx');
    fprintf("\nSeznam vyrovnaných měření:\nSměry:\nStanovisko   Cíl   směr[g]  Oprava[mgon]\n")
    fprintf("   %d       %d %8.4f   %8.5f\n",[smery(:,1:2),smery(:,3)-5,v(1:size(smery),1,:)/G2R*1000]')
    fprintf("Délky:\nStanovisko   Cíl   délka[m]  Oprava[mm]\n")
    fprintf("   %d       %d   %8.3f   %6.3f\n",[delky(:,1:2),delky(:,3),v(1:size(delky),1,:)*1000]')
else    
    fprintf("Nejsi v zenu bohužel /|\\ ^._.^ /|\\\n")
    fprintf("\nAposteriorní sm. odch. jednotková: %4.1f\n",s0);
    fprintf("%4.1f <= %4.1f <= %4.1f\n",[mez1,s0,mez2]);
    fprintf("\nSeznam vyrovnaných souřadnic:\n CB.    Y[m]         X[m]       OP[g]\n")
    fprintf("%d %12.3f %12.3f %8.4f\n",SS');
    fprintf("Seznam směrodatných odchylek vyrovnaných souřadnic:\n CB.  Y[mm]     X[mm]    OP[mgon]\n")
    fprintf("%d  %5.2f     %5.2f     %5.2f\n",Mx');
    fprintf("\nSeznam vyrovnaných měření:\nSměry:\nStanovisko   Cíl   směr[g]  Oprava[mgon]\n")
    fprintf("   %d       %d %8.4f   %8.5f\n",[smery(:,1:2),smery(:,3)-5,v(1:size(smery),1,:)/G2R*1000]')
    fprintf("Délky:\nStanovisko   Cíl   délka[m]  Oprava[mm]\n")
    fprintf("   %d       %d   %8.3f   %6.3f\n",[delky(:,1:2),delky(:,3),v(1:size(delky),1,:)*1000]')
    websave("NO.mp3","https://www.myinstants.com/media/sounds/nooo.mp3");
    [y, Fs] = audioread("NO.mp3");
    PO=audioplayer(y,Fs);
    playblocking(PO)
end
%% funkce
function [X] = Prevod_SS_M_na_V(SS)
X=[];
for n=1:size(SS,1)
    X=[X,SS(n,2:end)];
end
end
function [D] = Prevod_SS_V_na_M(X,SS)
    D=[SS(:,1),reshape(X,3,4)'];
end
function [lc,lp,lm] = vypocet_l(smery,delky,SS)
    G2R=pi/200;
    for n=1:size(smery,1)
        st=find(smery(n,1)==SS(:,1));
        cl=find(smery(n,2)==SS(:,1));
        lp(n)=atan2(SS(cl,2)-SS(st,2),SS(cl,3)-SS(st,3));
        if lp(n)<0
           lp(n)=lp(n)+2*pi;
        end
        lp(n)=lp(n)-SS(st,4);
        if lp(n)<0
           lp(n)=lp(n)+2*pi;
        end
    end
    m=size(lp,2);
    for n=1:size(delky,1)
        st=find(delky(n,1)==SS(:,1));
        cl=find(delky(n,2)==SS(:,1));
        lp(n+m)=sqrt((SS(cl,2)-SS(st,2))^2+((SS(cl,3)-SS(st,3)))^2);
    end
    lp=lp';
    lm=[smery(:,3)*G2R;delky(:,3)];
    lc=lm-lp;
end
function [B,b,b1] = Matice_B_b(zach,SS,X)
    B=zeros(3,length(X));
    b=zeros(3,3);
    st=find(zach(1,1)==SS(:,1));
    cl=find(zach(1,2)==SS(:,1));
    sm=atan2(SS(cl,2)-SS(st,2),SS(cl,3)-SS(st,3));
    if sm<0
       sm=sm+2*pi;
    end
    Y1=SS(st,2);
    X1=SS(st,3);
    PBYS=find(Y1==X);
    PBXS=find(X1==X);
    Y1=SS(cl,2);
    X1=SS(cl,3);
    PBYC=find(Y1==X);
    PBXC=find(X1==X);
    B(2,PBYS)=1;
    B(1,PBXS)=1;
    B(3,PBYC)=cos(sm);
    B(3,PBXC)=-sin(sm);
    b1=zeros(3,1);
end
function [A] = MAT_A(X,SS,smery,delky)
    A=zeros(size(smery,1)+size(delky,1),size(X,2));
    for n=1:size(smery,1)
        st=find(smery(n,1)==SS(:,1));
        cl=find(smery(n,2)==SS(:,1));
        Y1=SS(st,2);
        X1=SS(st,3);
        OP1=SS(st,4);
        AYST=find(Y1==X);
        AXST=find(X1==X);
        AOST=find(OP1==X);
        Y2=SS(cl,2);
        X2=SS(cl,3);
        AYCL=find(Y2==X);
        AXCL=find(X2==X);
        A(n,AYST)=1/(((Y1 - Y2)^2/(X1 - X2)^2 + 1)*(X1 - X2));
        A(n,AXST)=-(Y1 - Y2)/(((Y1 - Y2)^2/(X1 - X2)^2 + 1)*(X1 - X2)^2);
        A(n,AOST)=-1;
        A(n,AYCL)=-1/(((Y1 - Y2)^2/(X1 - X2)^2 + 1)*(X1 - X2));
        A(n,AXCL)=(Y1 - Y2)/(((Y1 - Y2)^2/(X1 - X2)^2 + 1)*(X1 - X2)^2);
    end
    for m=1:size(delky,1)
        st=find(delky(m,1)==SS(:,1));
        cl=find(delky(m,2)==SS(:,1));
        Y1=SS(st,2);
        X1=SS(st,3);
        AYST=find(Y1==X);
        AXST=find(X1==X);
        Y2=SS(cl,2);
        X2=SS(cl,3);
        AYCL=find(Y2==X);
        AXCL=find(X2==X);
        A(n+m,AYST)=(2*Y1 - 2*Y2)/(2*sqrt((X1 - X2)^2 + (Y1 - Y2)^2));
        A(n+m,AXST)=(2*X1 - 2*X2)/(2*sqrt((X1 - X2)^2 + (Y1 - Y2)^2));
        A(n+m,AYCL)=-(2*Y1 - 2*Y2)/(2*sqrt((X1 - X2)^2 + (Y1 - Y2)^2));
        A(n+m,AXCL)=-(2*X1 - 2*X2)/(2*sqrt((X1 - X2)^2 + (Y1 - Y2)^2));
    end
end


