A1=0.5612; 
A2=0.1403;
A3=0.5612;
x1=0;
x2=0.2113;
x3=1;
%convergent
Co=A1;
C2=3*(A2-A1)./(x2-x1).^2;
C3=-2*(A2-A1)./(x2-x1).^3;
X=linspace(0,0.2113,1000000);
A=Co+C2*(X-x1).^2+C3*(X-x1).^3;
D1=sqrt(A.*4./pi);

%plot(X,A)
%% divergent
do=A2;
d2=3.*(A3-A2)./(x3-x2).^2 ;
d3=-2.*(A3-A2)./(x3-x2).^3;
X2=linspace(0.2113,1,1000000);
A22=do+d2.*(X2-x2).^2+d3.*(X2-x2).^3;
%plot(X2,A2)
A2ls=-A22;
%A11=-A;
%XX=X +X2;
%Atotal=A+A22;
D2=sqrt(A22.*4./pi);
figure (1)
plot(X,D1,X,-D1,X2,D2,X2,-D2)
 %plot(XX,Atotal)
 grid minor 
 xlabel ('X axis')
 ylabel('diameter change')
Poin=3721.7.*10.^3;
Toin=840;
Min=0.137;
f=0.02;
Q=400.*10.^3;
gamma=1.4;

% convergent for mach change
dadx=2*C2*(X-x1)+3*C3*(X-x1).^2;
th=@(mt)[mt-(A1./A2)*(Min.*(1+(gamma-1).*Min.^2./2).^((-gamma-1)./(2*(gamma-1))))./(1+(gamma-1).*mt^2./2)^((-gamma-1)./(2*(gamma-1)))];
intialguess= 0;
MTh=fsolve(th,intialguess);
display(MTh)
M=linspace(0.137,MTh,1000);
%plot(X,M);
%% for convergent part
dx=(x2-x1)./999999; %spacing
da=dx*dadx;
i=1000;
R=287;
dq=Q.*dx./(1); %1m in length of duct
Cp=1005;
for i=1:999999;
    if i==1
        m(1)=Min;
        A11(1)=A1;
        Poini(1)=Poin;
        Po1(1)=Poin;
        Toini(1)=Toin;
        Tto(1)=Toin;
        m1f(1)=Min;
dm(1)=m(1).*((-2-(gamma-1).*m(1).^2)./(2.*(1-m(1).^2)))*da(1)./A11(1);

T(1)=Tto(1)./(1+(gamma-1).*m(1).^2./2);
dt(1)=-T(1).*((gamma-1).*m(1).^2)./((1+(gamma-1).*m(1)^2)./2).*dm(1)./(m(1));
P(1)=Po1(1).*(T(1)./Toini(1)).^(gamma./(gamma-1)); 
dp(1)=-P(1).*((gamma.*m(1).^2)./(1+((gamma-1).*m(1).^2)./2)).*dm(1)./m(1) ;
ROW(1)=P(1)./(R.*T(1)); %static density inlet 
ROWto(1)=Po1(1)./(Tto(1).*R); %total density inlet
a(1)=sqrt(gamma.*R.*T(1)); %speed of sound
v(1)=m(1).*a(1); %velocity at inlet
%for friction with area change
dm1(1)=m1f(1).*(gamma.*m1f(1).^2.*(1+((gamma-1).*m1f(1).^2)./2))./(2.*(1-m1f(1)^2)).*(f.*dx)./D1(1)+...
    m1f(1).*((-2-(gamma-1).*m1f(1).^2)./(2.*(1-m1f(1).^2)))*da(1)./A11(1);
Tto1f(1)=Toin; Po1f(1)=Poin;
T1f(1)=Tto1f(1)./(1+(gamma-1).*m1f(1).^2./2);
dtf(1)=-T1f(1).*(((gamma-1).*m1f(1).^2)./(1+((gamma-1).*m1f(1)^2)./2)).*dm1(1)./m1f(1);
P1f(1)=Po1f(1).*(T1f(1)./Tto1f(1)).^(gamma./(gamma-1));
dp1f(1)=-P1f(1).*((gamma.*m1f(1).^2)./(1+((gamma-1).*m1f(1).^2)./2)).*dm1(1)./m1f(1)...
    -P1f(1).*((gamma.*m1f(1).^2)./2).*(f.*dx)./D1(1);
Rowf(1)=P1f(1)./(R.*T1f(1));
Rowfo(1)=Po1f(1)./(R.*Tto1f(1));
af1(1)=sqrt(gamma.*R.*T1f(1));
vf1(1)=m1f(1).*af1(1);
% for friction & heat and area change
mt1(1)=Min; Tot(1)=Toin; Pot(1)=Poin;
dmt1(1)=mt1(1).*(((1+gamma.*mt1(1).^2).*(1+((gamma-1).*mt1(1).^2))./2)./(2.*(1-mt1(1).^2))).*(dq./(Cp.*Tot(1)))...
    +mt1(1).*(gamma.*mt1(1).^2.*(1+((gamma-1).*mt1(1).^2)./2))./(2.*(1-mt1(1)^2)).*(f.*dx)./D1(1)...
    +mt1(1).*((-2-(gamma-1).*mt1(1).^2)./(2.*(1-mt1(1).^2)))*da(1)./A11(1); %delta mach.
T1h(1)=Tot(1)./(1+((gamma-1).*mt1(1).^2)./2);%static temp.
dth(1)=-T1h(1).*(((gamma-1).*mt1(1).^2)./(1+((gamma-1).*mt1(1)^2)./2)).*dmt1(1)./mt1(1)...
    +(dq.*T1h(1))./(Cp.*Tot(1));% delta temp.
P1ft(1)=Pot(1).*(T1h(1)./Tot(1))^(gamma./(gamma-1));
dp1ft(1)=-P1ft(1).*((gamma.*mt1(1).^2)./(1+((gamma-1).*mt1(1).^2)./2)).*dmt1(1)./mt1(1)... %delta pressure
    -P1ft(1).*((gamma.*mt1(1).^2)./2).*(f.*dx)./D1(1)+...
    -P1ft(1).*(gamma.*mt1(1).^2.*dq)./(2.*Cp.*Tot(1));
rows(1)=P1ft(1)./(R.*T1h(1)); %static density 
rowso(1)=Pot(1)./(R.*Tot(1)); %tot. density
att(1)=sqrt(gamma.*R.*T1h(1)); %speed of sound
vv1(1)=att(1).*mt1(1); %velocity
    end
    dm(i+1)=m(i).*((-2-(gamma-1).*m(i).^2)./(2.*(1-m(i).^2)))*da(i)./A11(i);
    m(i+1)=dm(i)+m(i);
    A11(i+1)=da(i)+A11(i);
    dt(i+1)=-T(i).*(((gamma-1).*m(i).^2)./(1+((gamma-1).*m(i)^2)./2)).*dm(i)./m(i);
    T(i+1)=T(i)+dt(i);
    Tto(i+1)=T(i+1).*(1+((gamma-1).*m(i+1).^2)./2);
    dp(i+1)=-P(i).*((gamma.*m(i).^2)./(1+((gamma-1).*m(i).^2)./2)).*dm(i)./m(i) ;
    P(i+1)=P(i)+dp(i);
    Po1(i+1)=P(i+1).*(Tto(i+1)./T(i+1)).^(gamma./(gamma-1));
   ROW(i+1)=P(i+1)./(R.*T(i+1));
  ROWto(i+1)=Po1(i+1)./(Tto(i+1).*R);
   a(i+1)=sqrt(gamma.*R.*T(i+1));
   v(i+1)=m(i+1).*a(i+1);
   %for friction and area change
  dm1(i+1)=m1f(i).*(gamma.*m1f(i).^2.*(1+((gamma-1).*m1f(i).^2)./2))./(2.*(1-m1f(i)^2)).*(f.*dx)./D1(i)+...
       m1f(i).*((-2-(gamma-1).*m1f(i).^2)./(2.*(1-m1f(i).^2)))*da(i)./A11(i);
  m1f(i+1)=dm1(i)+m1f(i);
  dtf(i+1)=-T1f(i).*(((gamma-1).*m1f(i).^2)./(1+((gamma-1).*m1f(i)^2)./2)).*dm1(i)./m1f(i);
   T1f(i+1)=T1f(i)+dtf(i);
   Tto1f(i+1)=T1f(i+1).*(1+((gamma-1).*m1f(i+1).^2)./2);
   dp1f(i+1)=-P1f(i).*((gamma.*m1f(i).^2)./(1+((gamma-1).*m1f(i).^2)./2)).*dm1(i)./m1f(i)...
    -P1f(i).*((gamma.*m1f(i).^2)./2).*(f.*dx)./D1(i);
P1f(i+1)=P1f(i)+dp1f(i);
Po1f(i+1)=P1f(i+1).*(Tto1f(i+1)./T1f(i+1)).^(gamma./(gamma-1));
Rowf(i+1)=P1f(i+1)./(R.*T1f(i+1));
Rowfo(i+1)=Po1f(i+1)./(R.*Tto1f(i+1));
af1(i+1)=sqrt(gamma.*R.*T1f(i+1));
vf1(i+1)=m1f(i+1).*af1(i+1);
% for friction & area & heat
dmt1(i+1)=mt1(i).*(((1+gamma.*mt1(i).^2).*(1+((gamma-1).*mt1(i).^2)./2))./(2.*(1-mt1(i).^2))).*(dq./(Cp.*Tot(i)))...
    +mt1(i).*((gamma.*mt1(i).^2.*(1+((gamma-1).*mt1(i).^2)./2))./(2.*(1-mt1(i)^2))).*(f.*dx)./D1(i)...
    +mt1(i).*((-2-(gamma-1).*mt1(i).^2)./(2.*(1-mt1(i).^2)))*da(i)./A11(i);
mt1(i+1)=mt1(i)+dmt1(i);
dth(i+1)=-T1h(i).*(((gamma-1).*mt1(i).^2)./(1+((gamma-1).*mt1(i)^2)./2)).*dmt1(i)./mt1(i)...
    +(dq.*T1h(i))./(Cp.*Tot(i));
T1h(i+1)=T1h(i)+dth(i);
Tot(i+1)=T1h(i+1).*(1+((gamma-1).*mt1(i+1).^2)./2);
dp1ft(i+1)=-P1ft(i).*((gamma.*mt1(i).^2)./(1+((gamma-1).*mt1(i).^2)./2)).*dmt1(i)./mt1(i)...
    -P1ft(i).*((gamma.*mt1(i).^2)./2).*(f.*dx)./D1(i)+...
    -P1ft(i).*(gamma.*mt1(i).^2.*dq)./(2.*Cp.*Tot(i));
P1ft(i+1)=dp1ft(i)+P1ft(i);
Pot(i+1)=P1ft(i+1).*(1+((gamma-1).*mt1(i+1).^2)./2)^(gamma./(gamma-1));
rows(i+1)=P1ft(i+1)./(R.*T1h(i+1));
rowso(i+1)=Pot(i+1)./(R.*Tot(i+1));
att(i+1)=sqrt(gamma.*R.*T1h(i+1));
vv1(i+1)=att(i+1).*mt1(i+1);
end
%%
% figure (2)
 %plot(X,m)
mth=m(1000000) %mach at throat 
ath=A11(1000000) %area at throat 
Tth=T(1000000)  %static temp. at throat
pth=P(1000000)   %static pressure at throat
ROWth=ROW(1000000) %static density at throat
vth=v(1000000) %velocity at throat
Toth=Tto(1000000) % total temp. at throat
Poth=Po1(1000000) % total pressue at throat
m1fc=m1f(1000000)% mach at throat due to friction & area change
Tthf=T1f(1000000) %static temp. at throat due to friction& area change
Tothf=Tto1f(1000000) %total temp. at throat due to friction& area change
Pthf=P1f(1000000)
Pothf=Po1f(1000000)
mtth=mt1(1000000)
TShf=T1h(1000000)
toth=Tot(1000000)
poth=Pot(1000000)
ros=rows(1000000)
roso=rowso(1000000)
ato=att(1000000)
vvv=vv1(1000000)
m1=mth;
   %%  for divergent part
dadx2=2*d2*(X2-x2)+3*d3*(X2-x2).^2;
dx2=(x3-x2)./999999;
da2=dadx2*dx2;
dq2=Q.*dx2./(1);

for i=1:999999;
    if i==1
        mlu(1)=mth;
        Alu(1)=ath;
        Poin2(1)=Poth;
        Toin(1)=Toin;
        Tto2(1)=Toth;
dmlu(1)=mlu(1).*((-2-(gamma-1).*mlu(1).^2)./(2.*(1-mlu(1).^2)))*da2(1)./Alu(1);
T2(1)=Tto2(1)./(1+(gamma-1).*mlu(1).^2./2);
dt2(1)=-T2(1).*((gamma-1).*mlu(1).^2)./((1+(gamma-1).*mlu(1)^2)./2).*dmlu(1)./(mlu(1));
P2(1)=Poin2(1).*(T2(1)./Tto2(1)).^(gamma./(gamma-1)); %assume Po const as no NSW occur
dp2(1)=-P2(1).*((gamma.*mlu(1).^2)./(1+((gamma-1).*mlu(1).^2)./2)).*dmlu(1)./mlu(1) ;
ROW2(1)=P2(1)./(R.*T2(1));  %static density after throat 
ROWto2(1)=Poin2(1)./(Tto2(1).*R);
alu(1)=sqrt(gamma.*R.*T2(1)); %speed of sound
vlu(1)=mlu(1).*alu(1); %velocity at inlet
% due to friction &area change
 m1f2(1)=m1fc;
 T2f(1)=Tthf;
 Tto2f(1)=Tothf; Po2f(1)=Pothf;
dm2(1)=m1f2(1).*(gamma.*m1f2(1).^2.*(1+((gamma-1).*m1f2(1).^2)./2))./(2.*(1-m1f2(1)^2)).*(f.*dx2)./D2(1)+...
    m1f2(1).*((-2-(gamma-1).*m1f2(1).^2)./(2.*(1-m1f2(1).^2)))*da2(1)./Alu(1);
T2f(1)=Tto2f(1)./(1+(gamma-1).*m1f2(1).^2./2);
dtf2(1)=-T2f(1).*(((gamma-1).*m1f2(1).^2)./(1+((gamma-1).*m1f2(1)^2)./2)).*dm2(1)./m1f2(1);
P2f(1)=Po2f(1).*(T2f(1)./Tto2f(1)).^(gamma./(gamma-1));
dp2f(1)=-P2f(1).*((gamma.*m1f2(1).^2)./(1+((gamma-1).*m1f2(1).^2)./2)).*dm2(1)./m1f2(1)...
    -P2f(1).*((gamma.*m1f2(1).^2)./2).*(f.*dx2)./D2(1);
Rowf2(1)=P2f(1)./(R.*T2f(1));
Rowf2o(1)=Po2f(1)./(R.*Tto2f(1));
af2(1)=sqrt(gamma.*R.*T2f(1));
vf2(1)=m1f2(1).*af2(1);
% for friction & heat and area change
mt2(1)=mtth; Tot2(1)=toth; Pot2(1)=poth;
T2h(1)=TShf;
dmt2(1)=mt2(1).*(((1+gamma.*mt2(1).^2).*(1+((gamma-1).*mt2(1).^2))./2)./(2.*(1-mt2(1).^2))).*(dq2./(Cp.*Tot2(1)))...
    +mt2(1).*((gamma.*mt2(1).^2.*(1+((gamma-1).*mt2(1).^2)./2))./(2.*(1-mt2(1)^2))).*(f.*dx2)./D2(1)...
    +mt2(1).*((-2-(gamma-1).*mt2(1).^2)./(2.*(1-mt2(1).^2)))*da2(1)./Alu(1);

%T2h(1)=Tot2(1)./(1+(gamma-1).*mt2(1).^2./2);
dth2(1)=-T2h(1).*(((gamma-1).*mt2(1).^2)./(1+((gamma-1).*mt2(1)^2)./2)).*dmt2(1)./mt2(1)...
    +(dq2.*T2h(1))./(Cp.*Tot2(1));
P1ft2(1)=Pot2(1).*(T2h(1)./Tot2(1))^(gamma./(gamma-1));
dp1ft2(1)=-P1ft2(1).*((gamma.*mt2(1).^2)./(1+((gamma-1).*mt2(1).^2)./2)).*dmt2(1)./mt2(1)...
    -P1ft2(1).*((gamma.*mt2(1).^2)./2).*(f.*dx2)./D2(1)+...
    -P1ft2(1).*(gamma.*mt2(1).^2.*dq2)./(2.*Cp.*Tot2(1));
rowss(1)=ros; rowsoo(1)=roso;
atto(1)=ato; vv2(1)=vvv;
    end
    dmlu(i+1)=mlu(i).*((-2-(gamma-1).*mlu(i).^2)./(2.*(1-mlu(i).^2)))*da2(i)./Alu(i);
    mlu(i+1)=dmlu(i)+mlu(i);
    Alu(i+1)=da2(i)+Alu(i);
    dt2(i+1)=-T2(i).*(((gamma-1).*mlu(i).^2)./(1+((gamma-1).*mlu(i)^2)./2)).*dmlu(i)./mlu(i);
    T2(i+1)=T2(i)+dt2(i);
     Tto2(i+1)=T2(i).*(1+((gamma-1).*mlu(i).^2)./2);
     dp2(i+1)=-P2(i).*((gamma.*mlu(i).^2)./(1+((gamma-1).*mlu(i).^2)./2)).*dmlu(i)./mlu(i) ;
    P2(i+1)=P2(i)+dp2(i);
    Poin2(i+1)=P2(i+1).*(Tto2(i+1)./T2(i+1)).^(gamma./(gamma-1));
     ROW2(i+1)=P2(i+1)./(R.*T2(i+1));
     ROWto2(i+1)=Poin2(i+1)./(Tto2(i+1).*R);
     alu(i+1)=sqrt(gamma.*R.*T2(i+1));
   vlu(i+1)=mlu(i+1).*alu(i+1);
   %for friction with area change 
   dm2(i+1)=m1f2(i).*(gamma.*m1f2(i).^2.*(1+((gamma-1).*m1f2(i).^2)./2))./(2.*(1-m1f2(i)^2)).*(f.*dx2)./D2(i)+...
       m1f2(i).*((-2-(gamma-1).*m1f2(i).^2)./(2.*(1-m1f2(i).^2)))*da2(i)./Alu(i);
   m1f2(i+1)=m1f2(i)+dm2(i);
dtf2(i+1)=-T2f(i).*(((gamma-1).*m1f2(i).^2)./(1+((gamma-1).*m1f2(i)^2)./2)).*dm2(i)./m1f2(i);
 T2f(i+1)=dtf2(i)+T2f(i);
 Tto2f(i+1)=T2f(i+1).*(1+((gamma-1).*m1f2(i+1).^2)./2);
 dp2f(i+1)=-P2f(i).*((gamma.*m1f2(i).^2)./(1+((gamma-1).*m1f2(i).^2)./2)).*dm2(i)./m1f2(i)...
    -P2f(i).*((gamma.*m1f2(i).^2)./2).*(f.*dx2)./D2(i);
P2f(i+1)=P2f(i)+dp2f(i);
Po2f(i+1)=P2f(i+1).*(Tto2f(i+1)./T2f(i+1)).^(gamma./(gamma-1));
Rowf2(i+1)=P2f(i+1)./(R.*T2f(i+1));
Rowf2o(i+1)=Po2f(i+1)./(R.*Tto2f(i+1));
af2(i+1)=sqrt(gamma.*R.*T2f(i+1));
vf2(i+1)=m1f2(i+1).*af2(i+1);
% friction & area change &heat
dmt2(i)=mt2(i).*(((1+gamma.*mt2(i).^2).*(1+((gamma-1).*mt2(i).^2)./2))./(2.*(1-mt2(i).^2))).*(dq2./(Cp.*Tot2(i)))...
    +mt2(i).*((gamma.*mt2(i).^2.*(1+((gamma-1).*mt2(i).^2)./2))./(2.*(1-mt2(i)^2))).*(f.*dx2)./D2(i)...
    +mt2(i).*((-2-(gamma-1).*mt2(i).^2)./(2.*(1-mt2(i).^2)))*da2(i)./Alu(i);

mt2(i+1)=mt2(i)+dmt2(i);
dth2(i+1)=-T2h(i).*(((gamma-1).*mt2(i).^2)./(1+((gamma-1).*mt2(i)^2)./2)).*dmt2(i)./mt2(i)...
    +(dq2.*T2h(i))./(Cp.*Tot2(i));
T2h(i+1)=T2h(i)+dth2(i);
Tot2(i+1)=T2h(i+1).*(1+((gamma-1).*mt2(i+1).^2)./2);
dp1ft2(i+1)=-P1ft2(i).*((gamma.*mt2(i).^2)./(1+((gamma-1).*mt2(i).^2)./2)).*dmt2(i)./mt2(i)...
    -P1ft2(i).*((gamma.*mt2(i).^2)./2).*(f.*dx2)./D2(i)+...
    -P1ft2(i).*(gamma.*mt2(i).^2.*dq2)./(2.*Cp.*Tot2(i));
P1ft2(i+1)=P1ft2(i)+dp1ft2(i);
Pot2(i+1)=P1ft2(i+1).*(Tot2(i+1)./T2h(i)).^(gamma./(gamma-1));
rowss(i+1)=P1ft2(i+1)./(R.*T2h(i+1));
rowsoo(i+1)=Pot2(i+1)./(R.*Tot2(i+1));
atto(i+1)=sqrt(R.*gamma.*T2h(i+1));
vv2(i+1)=atto(i+1).*mt2(i+1);
end
%%
ae=Alu(1000000)
figure (3)
plot(X,m,X2,mlu,X,m1f,X2,m1f2,'-',X,mt1,X2,mt2)
xlabel('x-axis')
ylabel('mach no.')
grid minor 
legend ('area change conv.','area change div.','area & friction conv.','area &friction div.','A &F &Q conv.','A&F&Q div.')
figure (4)
plot(X,T,X2,T2,X,T1f,X2,T2f,X,T1h,X2,T2h)
grid minor 
xlabel('x-axis')
ylabel('Static T (Kelvin)')
legend ('area change conv.','area change div.','area & friction conv.','area &friction div.','A &F &Q conv.','A&F&Q div.','Location','Best')
figure (5)
plot(X,P,X2,P2,X,P1f,X2,P2f,X,P1ft,X2,P1ft2)
grid minor 
xlabel('x-axis')
ylabel('static P(pascal)')
legend ('area change conv.','area change div.','area & friction conv.','area &friction div.','A &F &Q conv.','A&F&Q div.','Location','Best')
figure (6)
plot(X,ROW,X2,ROW2,X,Rowf,X2,Rowf2,X,rows,X2,rowss)
grid minor 
xlabel('x-axis')
ylabel('static density (Kg./m^3)')
legend ('area change conv.','area change div.','area & friction conv.','area &friction div.','A &F &Q conv.','A&F&Q div.','Location','Best')
figure (7)
plot(X,ROWto,X2,ROWto2,X,Rowfo,X2,Rowf2o,X,rowso,X2,rowsoo)
xlabel('X-axis')
ylabel('total density (Kg./m^3)')
grid minor
legend ('area change conv.','area change div.','area & friction conv.','area &friction div.','A &F &Q conv.','A&F&Q div.','Location','Best')
figure (8)
plot(X,Tto,X2,Tto2,X,Tto1f,X2,Tto2f,X,Tot,X2,Tot2)
grid minor 
xlabel ('X-axis')
ylabel('total temp. (Kelvin)')
legend ('area change conv.','area change div.','area & friction conv.','area &friction div.','A &F &Q conv.','A&F&Q div.','Location','Best')
figure (9)
plot(X,Po1,X2,Poin2,X,Po1f,X2,Po2f,X,Pot,X2,Pot2)
grid minor
xlabel ('X-axis')
ylabel('total Pressure (pascal)')
legend ('area change conv.','area change div.','area & friction conv.','area &friction div.','A &F &Q conv.','A&F&Q div.','Location','Best')
figure (10)
plot(X,v,X2,vlu,X,vf1,X2,vf2,X,vv1,X2,vv2)
grid minor
xlabel('X-axis')
ylabel('velocity (m/s)')
legend ('area change conv.','area change div.','area & friction conv.','area &friction div.','A &F &Q conv.','A&F&Q div.','Location','Best')
figure (11)
plot(X,a,X2,alu,X,af1,X2,af2,X,att,X2,atto)
grid minor 
xlabel('X-axis')
ylabel('speed of sound (m/s)')
legend ('area change conv.','area change div.','area & friction conv.','area &friction div.','A &F &Q conv.','A&F&Q div.','Location','Best')
