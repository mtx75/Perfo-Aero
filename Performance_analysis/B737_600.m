clear all
clc
Vinf=100:1300;

Wmax=145498; %max take off weight
We=80031; %empty weight
Wav=(Wmax+We)./2;
Wmax2=125000 ; %max landing weight
rho0=0.002377; %density at sea level
rho1=8.91.*10.^-4; %density at 30,000ft above sea level
t0=518.67; %temperature at sea level
t1=411.84; %temperature at 30,000 ft 
b=112.53; %wing span 
S=1345.49; %wing platform area
Cl=(2.*Wav)./(rho1.*Vinf.^2.*S);
Cd0=0.02; %assumption 
AR=b.^2./S; %aspect ratio
e=0.77211; %oswald efficiency factor 
k3=1./(pi.*e.*AR);
k1=k3./3;
k=k1+k3;
k=k1+k3; %form factor
Cd= Cd0+k.*Cl.^2;
%%
figure (1)
  plot(Cd,Cl,'linewidth',1.5)
  xlabel('drag coeffiecient')
  ylabel('lift coeffiecient')
 grid on
  title('drag polar')
  axis([0 0.25 0 2])
C0=Cl./Cd; %lift to drag coeffecient
C1=(Cl.^0.5)./(Cd);
C2=(Cl.^1.5)./(Cd);
c0=560*ones(length(0:15));
c1=740*ones(length(0:22));
c2=425*ones(length(0:13));
%%
figure(2)
   
   plot(Vinf,C0,Vinf,C1,Vinf,C2 ,c0,0:15,c1,0:22,c2,0:13,'linewidth',1.5)
   title('V with ratio of Coeffecient of lift with drag')
   xlabel('V'); 
   ylabel('CL/CD,CL.^(0.5)/CD,CL.^(1.5)/CD')
     grid on
   legend('Cl/Cd','Cl^0^.^5/Cd','Cl^1^.^5/Cd','Cl^1^.^5/Cd max','Cl/Cd max','Cl^0^.^5/Cd max');
Tr=0.5.*rho1.*S.*Vinf.^2.*Cd; %thrust required
r=1716; %universal gas constant
y=1.4;
N=2;   %no of engines
Ao=0.2; %assumption
n=0.6; %assumption
Tavo=24200; %static thrust
a1=sqrt(y.*r.*t1); %speed of sound @30000 ft
Tav1=Tavo.*(rho1./rho0).^n;     
Ta=Tavo.*N.*Ao.*Vinf.^-0.6.*a1.^0.6;
Vstall=sqrt((2.*Wav)./(rho1.*S.*2.46));
vstall=Vstall*ones(1,length(0:1.718*10^4));
Vmax=861; 
V2nd=236;  
v2nd=V2nd*ones(1,length(0:2.291*10^4));
vmax=Vmax*ones(1,length(0:1.056*10^4));

%%
figure (3)

   plot(Vinf,Tr,Vinf,Ta,vstall,0:1.718*10^4,vmax,0:1.056*10^4,v2nd,0:2.291*10^4,'linewidth',1.5)
   xlabel('velocity (ft/sec)')
   ylabel('thrust (lb)')
   title('thrust vs velocity')
   grid on
   legend('Thrust Required','Thrust Available','Vstall','Vmax','V2nd');

   
 Pr=(Tr.*Vinf)./550;
 Pa=(Ta.*Vinf)./550;
 vstall1=Vstall*ones(1,length(0:0.8620*10^4));
 v2nd1=V2nd*ones(1,length(0:9832));
 vmax1=Vmax*ones(1,length(0:1.653*10^4));
 
 figure (4)
     plot(Vinf,Pr,Vinf,Pa,vstall1,0:0.8620*10^4,vmax1,0:1.653*10^4,v2nd1,0:9832,'linewidth',1.5)
     xlabel('velocity ft/sec')
     ylabel('power hp')
     title('velocity vs power')
       grid on
     legend('Power Required','Power Available','Vstall','Vmax','V2nd point');
Vinf2=350:1300;     
Cl2=(2.*Wav)./(rho1.*Vinf2.^2*S);     
Cd2= Cd0+k.*Cl2.^2;
Tr2=0.5.*rho1.*S.*Vinf2.^2.*Cd2;
Ta2=Tavo.*N.*Ao.*Vinf2.^-0.6.*a1.^0.6;
     
B=(Ta2-Tr2)./Wav;
C=Vinf2.*B;
Theta=asin(C./Vinf2);
D=Vinf2.*cosd(Theta);

Vinf4=0:500;
Y=0.06296.*Vinf4;
Vinf5=0:533;
Ratemax=0.0559475.*Vinf5;
%%
figure (5)
    plot(D,C,Vinf4,Y,'--',Vinf5,Ratemax,'linewidth',1.5)
    
    xlabel('horizontal velocity ft/sec')
    ylabel('vertical velocity ft/sec')
     grid on
    title('Vh vs Vv')
    legend('Hodograph for Climbing','V theta Maximum','Rate of Climb Maximum');
    axis([0 1200 0 80])
    %%
B2=(Ta-Tr)./Wav;
C2=Vinf.*B2;
Theta=asin(B2);
   
theta1=(Theta.*180)./pi;
Vmaxtheta=348;
vmaxtheta=Vmaxtheta*ones(length(0:11));

Vv1=(Pr.*550)./Wav; %assume weight tends to empty weight of the A/C

figure (6)
   plot(Vinf,Vv1,'linewidth',1.5)
   Vv1=gca ;
   Vv1.XAxisLocation='top' ;
   Vv1.YDir='reverse' ;
   xlabel('velocity free stream ')
   ylabel('Rate of decent')
   title('Rate of decent vs free stream velocity')
     grid on
    
Vv2=(Pr.*550)./Wav;
F=Vv2./Vinf;
theta2=asin(F);
Vh1=Vinf.*cosd(theta2);
Vh2=0:560;
m1=0.068325;
Z=Vh2.*m1;
Y=33.5*ones(1,length(0:450));
Vh3=0:440;
m2=0.0763098;
Z2=m2.*Vh3;

%%
figure (7)
     plot(Vh1,Vv2,Vh2,Z,0:450,Y,'--',Vh3,Z2,'linewidth',2)
    Vv2=gca ;
    Vv2.XAxisLocation='top' ;
    Vv2.YDir='reverse' ;
    xlabel('horizontal velocity')
    ylabel('rate of decent')
    title('Hodograph for Gliding')
    legend('Hodograph for gliding','V(theta) Minimum','V(v) Minimum','V(r/d) Minimum');
      grid on
h= 0:43000;    
Z=42393.6098./-836.1659;
U=1./-836.1659;
h2=30000;

Rate_climb=-Z+U.*h;
Rate_service=100./60*ones(1,length(0:4.1*10^4));
figure(8)
    plot(Rate_climb,h,Rate_service,0:4.1*10^4,'--' ,'linewidth',1.7);
    ylabel('Altitude')
    xlabel('Rate of climb')
    title ( 'Rate of climb VS Altitude')
    axis([0 50 0 45000])
      grid on
    legend('Rate of climb VS Altitude','Service ceiling');
tmin=(log(-Z+U.*h2)-log(-Z))./U;
%%
 figure (9)
  Vinf3=0:750;
    m=11.924;
    T=Vinf3.*m;
    X=7482*ones(1,length(300:620));
   plot(Vinf,Tr,Vinf3,T,300:620,X,'linewidth',2)
   legend('Thrust Required','Range Max','Endurance Max');
  xlabel('Velocity')
  ylabel('Thrust Required')
  title('Range and Endurance')
    grid on
    %%
Ct=0.4/3600; %pound/hr
      Cl__Cd_max=14.63;
rho1=8.91.*10.^-4; %density at 30,000ft above sea level

Range_max1=2./Ct.*sqrt(2./rho1./S).*21.8.*(Wmax.^0.5-Wmax2.^0.5) ; % ft %at h=30000
E_max1=1./Ct.*Cl__Cd_max.*log(Wmax./Wmax2) ; % seconds %at h=30000
rho_max=5.87.*10.^-4; %density at 40,000ft above sea level   
Range_max2=2./Ct.*sqrt(2./rho_max./S).*21.8.*(Wmax.^0.5-Wmax2.^0.5) ; % ft %at h=40000
E_max2=1./Ct.*Cl__Cd_max.*log(Wmax./Wmax2) ; % seconds?%at?h=40000
%% level turn
a0=sqrt(y.*r.*t0);
c=Tavo.*N.*0.4.*a0.^0.6;
g=32.17;
%% Maximum maximum load factor
Vn_M=[((2-n).*(c./Wmax).*(Wmax./S))./(2.*rho0.*Cd0)].^(1./(n+2))
n_M=((0.5*rho0.*Vn_M^2)./(k*(Wmax/S))*(c.*Vn_M^-0.6/Wmax - 0.5.*rho0.*Vn_M^2*Cd0/(Wmax./S)))^0.5
Rn_M=Vn_M.^2./(g.*sqrt(n_M.^2-1))
omegan_M=(g.*sqrt(n_M.^2-1))./Vn_M
CLmax=2.46;

Vrmin=((4.*k*(Wmax/S))./(rho0.*(c./Wmax).*(1+(n/2)))).^(1./(2-n))
n_vrmin=((0.5*rho0.*Vrmin^2)./(k*(Wmax/S))*(c.*Vrmin^-0.6/Wmax - 0.5.*rho0.*Vrmin^2*Cd0/(Wmax./S)))^0.5
CLrmin=2*n_vrmin*Wmax./(rho0*Vrmin.^2.*S)

%% R min
stall=@(ViStall)[ ViStall.^2 - (2.*c.*ViStall.^-0.6./S)./(rho0.*(Cd0+k.*CLmax.^2))];
    intialguess= 100;
    VinfStall= fsolve(stall,intialguess);
    disp(VinfStall)
    
if CLrmin<CLmax
Rn_vrmin=Vrmin.^2./(g.*sqrt(n_vrmin.^2-1))
omega_rmin=(g.*sqrt(n_vrmin.^2-1))./Vrmin
elseif CLrmin> CLmax
    Vrmin=VinfStall
    n_vrmin=((0.5*rho0.*Vrmin^2)./(k*(Wmax/S))*(c.*Vrmin^-0.6/Wmax - 0.5.*rho0.*Vrmin^2*Cd0/(Wmax./S)))^0.5
    Rn_vrmin=Vrmin.^2./(g.*sqrt(n_vrmin.^2-1))
omega_rmin=(g.*sqrt(n_vrmin.^2-1))./Vrmin
end
%% Omega Max

vomax=@(VOMax)[ VOMax-((4.*k*(Wmax/S).^2)./(rho0.^2.*Cd0)- ((c./S).*0.6.*VOMax.^(2-0.6)./(rho0.*Cd0))).^0.25];
    intialguess= 100;
    VOmax= fsolve(vomax,intialguess);
    disp(VOmax)
 %%
n_Omax=((0.5*rho0.*VOmax^2)./(k*(Wmax/S))*(c.*VOmax^-0.6/Wmax - 0.5.*rho0.*VOmax^2*Cd0/(Wmax./S)))^0.5
CLO=2*n_Omax*Wmax./(rho0*VOmax.^2.*S)
%%
if CLO<CLmax
Rn_vOmax=VOmax.^2./(g.*sqrt(n_Omax.^2-1))
omegaMax=(g.*sqrt(n_M.^2-1))./VOmax
elseif CLO> CLmax
    VOmax=VinfStall
    n_Omax=((0.5*rho0.*VOmax^2)./(k*(Wmax/S))*(c.*VOmax^-0.6/Wmax - 0.5.*rho0.*VOmax^2*Cd0/(Wmax./S)))^0.5
    Rn_vOmax=VOmax.^2./(g.*sqrt(n_Omax.^2-1))
omegaMax=(g.*sqrt(n_Omax.^2-1))./VOmax
end
%% load factor with velocity
vnn= 273:846;
n_max=((0.5.*rho0.*vnn.^2)./(k*(Wmax./S)).*(c.*vnn.^-0.6./Wmax - 0.5.*rho0.*vnn.^2*Cd0/(Wmax./S))).^0.5 ;
Stall=0:276;
n_maxstall=0.5.*rho0.*Stall.^2.*CLmax./(Wmax./S);
    figure (10)
plot(vnn,n_max,Stall, n_maxstall ,'linewidth',1.5)
grid minor
xlabel('Vinf ft/sec')
ylabel('n_max')
legend('n_max thrust','n_max stall')
%% load factor
vnnn=0:322;
n_maxstall2=0.5.*rho0.*vnnn.^2.*CLmax./(Wmax./S);
voper=2.79.*ones(length(322:846));
clmaxneg=-2; %assumption
nmaxneg=-1;
vstallneg=sqrt((-2.*Wmax)./(S.*rho0.*clmaxneg));
vneg=0:213;
n_maxstallneg2=0.5.*rho0.*vneg.^2.*clmaxneg./(Wmax./S)
nmaxneg=-1.*ones(length(213:846));
vneg2=846.*ones(length(-1:3));
figure (11)
plot (vnnn,n_maxstall2,322:846,voper,vneg,n_maxstallneg2,213:846,nmaxneg,vneg2,-1:3 ,'linewidth',1.5)
grid minor
xlabel('Vinf')
ylabel('nmax')
title('V_N')
axis([0 950 -2 4])
%% Takeoff Performance
rho0= 0.002377;

S= 1345.49;
Cl_maxTakeoff= 1.6.*cos(25);
%%
g= 32.17;
h_OB= 35;
Kuc= 4.485.*10.^-5;
Cd_0= 0.02;
h2 = 41;
b = 112.53;
Cl= 0.1;
mu1 = 0.04; %dry concrete brakes off

n= 0.6;
N= 3;

V_stall= sqrt((2./rho0).*(Wmax./S).*(1./Cl_maxTakeoff));
V_lo= 1.1.*V_stall;
V= 0.7.*V_lo;
R= (6.95.*(V_stall.^2))./g;
theta_OB= acos(1-(h_OB./R));

S_a_T= R.*sin(theta_OB);
delta_Cd_0= ((Wmax./S).*47.8778).*Kuc.*(Wmax.*0.4536).^-0.215;
G= ((16.*h2./b).^2)./(1+(16.*h2./b).^2);
K_A= -(rho0.*S)./(2.*Wmax).*(Cd_0 + delta_Cd_0 + (k1+ G*k3).*Cl.^2 - mu1.*Cl);

treal=51529.78857-38.859.*V+0.01988*V^2; %thrust for takeoff
K_T= (treal)./Wmax - mu1;
S_g_T=(1./(2.*g.*K_A)).*log(1+(K_A.*(V_lo.^2))./K_T)+ 3.*V_lo ;
S_tot_T= S_a_T + S_g_T;


%% landing
sweep=25;
Cl_max=cos(sweep).*(1.8+2.2)./2;
theta_a= 3; 

mu2 = 0.4; %brakes on dry concrete assumption
n= 0.6;
N= 3;
J= 1.15;
V_stall2= sqrt((2./rho0).*(Wmax2./S).*(1./Cl_max));

V_TD= 1.15.*V_stall2;
RL= [(1.23.*V_stall2).^2]./[0.2.*g];

S_f_L= RL.*sind(theta_a);
h_f= RL.*(1-cosd(theta_a));

S_a_L= (50-h_f)./tand(theta_a);
S_g_L=V_TD.*N+ (J.^2.*Wmax2)./(S.*g.*rho0.*Cl_max.*mu2);
S_tot_L= S_f_L + S_a_L + S_g_L;
%% landing 2
sweep=25;
Cl_max=cos(sweep).*(1.8+2.2)./2;
c=Tavo.*N.*0.4.*a0.^0.6;
theta_a= 2; 
KucL= 3.16 * 10^-5;
delta_Cd_0L= ((Wmax2./S).*47.8778).*KucL.*(Wmax2.*0.4536).^-0.215;
mu2 = 0.3; %brakes on dry concrete assumption range 0.3-0.5
n= 0.6;
N= 3;
J= 1.15;
V_stall2= sqrt((2./rho0).*(Wmax2./S).*(1./Cl_max));

V_TD= 1.15.*V_stall2; %range from 219 to 244
RL= [(1.23.*V_stall2).^2]./[0.2.*g];

S_f_L= RL.*sind(theta_a);
h_f= RL.*(1-cosd(theta_a));

S_a_L= (50-h_f)./tand(theta_a);
tlanding=0.3*c.*(0.7*V_TD).^-0.6; %30 percentage thrust reverse range 30-45
J_T= (tlanding)./(Wmax2) + mu2;
J_A= (rho0.*S)./(2.*Wmax2).*(Cd_0 + delta_Cd_0L + (k1+ G*k3).*(Cl).^2 - mu2.*Cl);
S_g_L2=(1./(2.*g.*J_A)).*log(1+(J_A.*(V_TD.^2))./J_T)+ 3.*V_TD ;

S_tot_L2= S_f_L + S_a_L + S_g_L2;
