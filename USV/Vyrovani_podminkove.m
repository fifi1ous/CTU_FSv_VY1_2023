clc; clear; format long g

% Měřená převýšení:

%     c.h.    h/m        d/m
zadani=[65   31.8140     326.3115    
       64   23.46216     750.3989
       62   -46.0387    1064.68565
       54   -8.31098    627.647
       24   69.54468    666.9285
       52   -77.839501    1136.005];

 %Přesnost měření [m]: 
  sh0 =  0.0044;  

 %Výchozí body [m]: 
  H106 =  [873.7409];  
  k=3;
  n=6;
  
  %%
  % Podminky
  syms h65 h64 h62 h54 h24 h52;
  
  prom = [h65 h64 h62 h54 h24 h52];
  
 eq1= h65 + h52 - h62;
 eq2= h54 - h24 - h52;
 eq3= h65 + h54 - h64;

 eqs = [eq1 eq2 eq3];
    
  %u=zeros(k,1);
  for i = 1:k
      u(i)=eval(subs(eqs(i),[prom],[zadani(:,2)']));
  end


  %%
  % Vahy
  
  sh=(zadani(:,3)/1000).^2;
  P=diag(1./sh);
  
  %%
  % Matice A
  for i=1:k
      for j=1:n
         A(j,i)=diff(eqs(i),prom(j));
      end
  end
  clear vars i j
  A=eval(A);
  
  
  
  
  %%
  % Výpočet normálních rovnic
  u=u';
  l=zadani(:,2);
  kol=-(inv(A'*inv(P)*A))*u;

  % Výpočet oprav
  v=inv(P)*A*kol;
  l_=l+v;

  s0=sqrt((v'*P*v)/(n-k));

%   Ml=s0^2*(inv(P)-inv(P)*A*(inv(A'*inv(P)*A))*A'*inv(P));
%   Ml_=sqrt(diag(Ml));

  Qh=s0^2*(A*inv(A'*P*A)*A');
  qh=sqrt(diag(Qh))';

  Qx=s0^2*(inv(A'*P*A));
  qx=sqrt(diag(Qx))';

  %% Výpočet výšek
  H105 = H106 + l_(1);
  H104 = H106 + l_(2);
  H102 = H106 + l_(3);

%   %%
%   % Výpočet nadmořských výšek
% ATH = [1 0 0 0 0 0
%        0 1 0 0 0 0
%        0 0 1 0 0 0];
%   
%   vysky2=H106+ATH*l_;
%   MH=s0^2*ATH*Ml*ATH';
%   MH_=sqrt(diag(MH));
% Aposteriorní směrodatná odchylka jednotková

kontrola=zeros(k,1);

for i=1:k
    kontrola(i)=eval(subs(eqs(i),[prom],[l_']));
end

%% Kvantily
alfa = 0.05;
kv1 = chi2inv(alfa/2,n-k);
kv2 = chi2inv(1-alfa/2,n-k);

sp = sqrt(kv1/(n-k)); s0s0_ = s0/sh0; hor = sqrt(kv2/(n-k));



%%
for i=1:n
    fprintf('Převýšení %d je %.5fm\n',zadani(i,1),l_(i))
    fprintf('Směrodatná odchylka převýšení %d je %.4fm\n\n',zadani(i,1),qh(i))
end

body=[105 104 102];
vysky=[H105,H104,H102];
for i=1:k
    fprintf('Vyska bodu %d je %.5fm\n',body(i),vysky(i))
    fprintf('Směrodatná odchylka vyrovnané výšky %d je %.4f\n\n',body(i),qx(i))
end

fprintf('Aposteriorní směrodatná odchylka jednotková je %.2f\n',s0)






  
  