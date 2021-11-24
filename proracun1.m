clc
clear all
close all

%%Pocetna pozicija

O1 = [0,0];
A1 = [-2732,-727];
C1 = [-550,-980];

O2 = [-6259,-4311];
A2 = [-6473,-3476];
C2 = [-4955,-2153];

H = [-5591,-11357];

G1 = [-3129; -2155.2];
G2 = [-5925;-7862];

GC1 = [];
GC2 = [];

% figure(10)
% axis ([-15000 9000 -13500 5500])
% hold on
% grid on
% plot(O1(1),O1(2),'k^')
% plot(A1(1),A1(2),'kx')
% plot(C1(1),C1(2),'k^')
% plot(O2(1),O2(2),'kx')
% plot(A2(1),A2(2),'kx')
% plot(C2(1),C2(2),'kx')
% line([O1(1), O2(1)], [O1(2), O2(2)], 'LineWidth',2)
% line([O2(1), H(1)], [O2(2), H(2)], 'LineWidth',2)
% plot(H(1),H(2),'kx')
% plot(G1(1),G1(2),'ro')
% plot(G2(1),G2(2),'ro')

ro11 = [O2(1)-O1(1);O2(2)-O1(2);0;];
ro22 = [H(1)-O2(1);H(2)-O2(2);0];

l1 = sqrt(ro11(1)^2+ro11(2)^2);
l2 = sqrt(ro22(1)^2+ro22(2)^2);
%% Pozicija tachaka

O1A1 = [A1(1)-O1(1);A1(2)-O1(2);0;];
O2A2 = [A2(1)-O2(1);A2(2)-O2(2);0];
O1C2 = [C2(1)-O1(1);C2(2)-O1(2);0];

O1GC1 = [-1641;-853;0];

O1G1 = [G1(1)-O1(1);G1(2)-O1(2);0];
O2G2 = [G2(1)-O2(1);G2(2)-O2(2);0];

la1 = sqrt(O1A1(1)^2+O1A1(2)^2);
la2 = sqrt(O2A2(1)^2+O2A2(2)^2);
lc2 = sqrt(O1C2(1)^2+O1C2(2)^2);

lg1 = sqrt(O1G1(1)^2+O1G1(2)^2);
lg2 = sqrt(O2G2(1)^2+O2G2(2)^2);

lgc1 = sqrt(O1GC1(1)^2+O1GC1(2)^2);
%% Inverzna kinematika

xH = -5591:50:-2871;
yH = 3.786*xH + 9811;

x1 = zeros;
y1 = zeros;
x2 = zeros;
y2 = zeros;
xa1 = zeros;
ya1 = zeros;
xa2 = zeros;
ya2 = zeros;
xc2 = zeros;
yc2 = zeros;
q1 = zeros;
q2 = zeros;

qmax = [76*2*pi/360; 90*2*pi/360];

q = [0;0];

for i = 1:length(xH)
    
    q2(i) = pi - acos((l1^2 + l2^2 - xH(i)^2 - yH(i)^2)/(2*l1*l2));
    q1(i) = atan2(yH(i),xH(i)) - atan2(l2*sin(q2(i)),(l1 + l2*cos(q2(i)))); 
  
x1(i) = l1*cos(q1(i));
y1(i) = l1*sin(q1(i));
x2(i) = x1(i) + l2*cos(q1(i) + q2(i));
y2(i) = y1(i) + l2*sin(q1(i) + q2(i));

xa1(i) = la1*cos(q1(i) - 0.3426);
ya1(i) = la1*sin(q1(i) - 0.3426);
%0.3426 ugao izmaedju O1A1 i O1O2

xc2(i) = lc2*cos(q1(i) - 0.1932);
yc2(i) = lc2*sin(q1(i) - 0.1932);
%0.1932 ugao izmaedju O1C2 i O1O2

xa2(i) = x1(i) + la2*cos(q1(i) + q2(i) - 2.9852);
ya2(i) = y1(i) + la2*sin(q1(i) + q2(i) - 2.9852);
%-2.9852 ugao izmaedju O2H i O2A2

xg1(i) = lg1*cos(q1(i));
yg1(i) = lg1*sin(q1(i));

xg2(i) = x1(i) + lg2*cos(q1(i) + q2(i));
yg2(i) = y1(i) + lg2*sin(q1(i) + q2(i));

xgc1(i) = (xa1(i) - C1(1))/2 + C1(1);
ygc1(i) = (ya1(i) - C1(2))/2 + C1(2);

xgc2(i) = (xa2(i) - xc2(i))/2 + xc2(i);
ygc2(i) = (ya2(i) - yc2(i))/2 + yc2(i);

end

figure(1)
hold on
axis ([-15000 9000 -13500 5500])
line([0, x1(1)], [0, y1(1)], 'LineWidth',2)
line([x1(1), x2(1)], [y1(1), y2(1)], 'LineWidth',2)
plot(0, 0, 'k^', 'MarkerSize',10, 'LineWidth',1.5)
plot(x1(1), y1(1), 'kx', 'MarkerSize',10, 'LineWidth',1.5)
plot(x2(1), y2(1), 'rx', 'MarkerSize',10, 'LineWidth',1.5)
plot(xH,yH,'k*')

plot(C1(1), C1(2), 'k^', 'MarkerSize',8, 'LineWidth',1.5)

xlabel('$$x$$','interpreter','latex','fontsize',14)
ylabel('$$y$$','interpreter','latex','fontsize',14)
%% Animacija
h = animatedline('Color','b','LineStyle','--','LineWidth',1.5); 
h2 = animatedline('LineStyle','none','Marker','o'); 
h3 = animatedline('LineStyle','none','Marker','o','MarkerSize',1.5,'Color','k');
h4 = animatedline('Color','r','LineWidth',1.5);
h5 = animatedline('Color','y','LineWidth',1.5);
h6 = animatedline('Color','r','LineWidth',1.5);
h7 = animatedline('Color','y','LineWidth',1.5);
h8 = animatedline('Color','g','Marker','o');
h9 = animatedline('Color','g','Marker','o');
h10 = animatedline('Color','g','Marker','o');
h11 = animatedline('Color','g','Marker','o');

n = 0.001;
for k = 1:length(xH)
    pause(n)
    addpoints(h,[0, x1(k)], [0, y1(k)])
    addpoints(h,[x1(k), x2(k)], [y1(k), y2(k)])
    addpoints(h2,xH(k),yH(k))
    
    addpoints(h4,[0, xa1(k)],[0, ya1(k)])
    addpoints(h6,[x1(k), xa2(k)],[y1(k), ya2(k)])
    
    addpoints(h5,[C1(1), xa1(k)], [C1(2), ya1(k)])
    addpoints(h7,[xc2(k), xa2(k)], [yc2(k), ya2(k)])
    
    addpoints(h,[x1(k), xc2(k)], [y1(k), yc2(k)])
    addpoints(h,[xc2(k), 0], [yc2(k), 0])
    
    addpoints(h8,xg1(k),yg1(k))
    addpoints(h9,xg2(k),yg2(k))
    
    addpoints(h10,xgc1(k),ygc1(k))
    addpoints(h11,xgc2(k),ygc2(k))
    
    drawnow
    pause(n)
    clearpoints(h)
    clearpoints(h2)
    clearpoints(h4)
    clearpoints(h5)
    clearpoints(h6)
    clearpoints(h7)
    clearpoints(h8)
    clearpoints(h9)
    clearpoints(h10)
    clearpoints(h11)
end
    line([0, x1(end)], [0, y1(end)], 'LineWidth',2)
    line([x1(end), x2(end)], [y1(end), y2(end)], 'LineWidth',2)
    plot(0, 0, 'kx', 'MarkerSize',10, 'LineWidth',1.5)
    plot(x1(end), y1(end), 'kx', 'MarkerSize',10, 'LineWidth',1.5)
    plot(x2(end), y2(end), 'rx', 'MarkerSize',10, 'LineWidth',1.5)  

figure(2)
subplot(2,2,1)
hold on
grid on
title('Zavisnost $$\varphi_2$$ od $$\varphi_1$$','interpreter','latex')
plot(rad2deg(q1),rad2deg(q2))
xlabel('$$\varphi_1 [^{\circ}]$$','interpreter','latex')
ylabel('$$\varphi_2 [^{\circ}]$$','interpreter','latex')

subplot(2,2,2)
hold on
grid on
title('Zavisnost  $$\Delta \varphi_2$$ od $$\Delta \varphi_1$$','interpreter','latex')
plot(rad2deg(abs(q1-q1(1))),rad2deg(q2-q2(1)))
xlabel('$$\Delta \varphi_1[^{\circ}]$$','interpreter','latex')
ylabel('$$\Delta \varphi_2[^{\circ}]$$','interpreter','latex')
%U opticaju vrednosti od -40 > 30, 153 > 63
%% Racunanje Sile

rB = zeros;
rBD = zeros;
rED = zeros;
rER = zeros;
rC = zeros;
rC2 = zeros;
lc1 = zeros;
lb = zeros;
lbd = zeros;
le = zeros; 
lER = zeros;
beta = zeros;
k1ER = zeros;
phi = zeros;
F = zeros;

for k = 1:length(xH) 
    
rB(k,1) = xa2(k) - x1(k);
rB(k,2) = ya2(k) - y1(k);       %A2 - O2

rBD(k,1) = xc2(k) - xa2(k);     %A2 - C2
rBD(k,2) = yc2(k) - ya2(k);   

rED(k,1) = x1(k) - xc2(k);      %A2 - C2
rED(k,2) = y1(k) - yc2(k);

rER(k,1) = xH(k) - x1(k); 
rER(k,2) = yH(k) - y1(k);

rC(k,1) = xa1(k) - C1(1);
rC(k,2) = ya1(k) - C1(2); 

rC2(k,1) = xc2(k) - xa2(k); %%%
rC2(k,2) = yc2(k) - ya2(k); 

lc1(k) = sqrt(rC(k,1).^2 + rC(k,2).^2);
lc2(k) = sqrt(rC2(k,1).^2 + rC2(k,2).^2);
lb(k) = sqrt(rB(k,1).^2 + rB(k,2).^2);
lbd(k) = sqrt(rBD(k,1).^2 + rBD(k,2).^2);
le(k) = sqrt(rED(k,1).^2 + rED(k,2).^2);
lER(k) = sqrt(rER(k,1).^2 + rER(k,2).^2);

beta(k) = acos((lb(k).^2 + lbd(k).^2 - le(k).^2)/(2.*lb(k).*lbd(k)));

k1ER(k) = (yH(k) - y1(k))/(xH(k) - x1(k));

k2 = 3.786;

phi(k) = atan(abs((k2-k1ER(k))/(1 + k1ER(k).*k2)));

FR = 3000;
%%%%%%%%%%%%%%%%%%%%
F(k) = lER(k)*1e-3*FR*sin(pi/2 + phi(k))/(lb(k)*1e-3*sin(pi - beta(k)));
%%%%%%%%%%%%%%%%%%%s
end

subplot(2,2,3)
hold on
grid on
title('Zavisnost sile od $$\Delta \varphi2$$','interpreter','latex')
plot(rad2deg(q2-q2(1)),F)
xlabel('$$\varphi_2 [^{\circ}]$$','interpreter','latex')
ylabel('$$F_R [N]$$','interpreter','latex')
%% Vrednost pritiska

dk = 140;
Ak = dk^2/4*pi*1e-6;
pk = F./Ak;
pk_b = pk./1e5;

subplot(2,2,4)
hold on
grid on
title('Zavisnost pritiska na strani klipa od $$\Delta \varphi2$$','interpreter','latex')
plot(rad2deg(q2-q2(1)),pk_b)
xlabel('$$\varphi_2 [^{\circ}]$$','interpreter','latex')
ylabel('$$p_k [bar]$$','interpreter','latex')
%% Izduezenje cilindara

figure(3)
subplot(1,2,1)
hold on
grid on
plot(rad2deg(abs(q1-q1(1))), lc1-lc1(1))
title('Izduzenje cilindra 1 od $$\Delta \varphi_1$$','interpreter','latex')
xlabel('$$\Delta \varphi_1 [^{\circ}]$$','interpreter','latex')
ylabel('$$l_1 [mm]$$','interpreter','latex')

subplot(1,2,2)
hold on
grid on
plot(rad2deg(q2-q2(1)), lc2-lc2(1))
title('Izduzenje cilindra 2 od $$\Delta \varphi_2$$','interpreter','latex')
xlabel('$$\Delta \varphi_2 [^{\circ}]$$','interpreter','latex')
ylabel('$$l_2 [mm]$$','interpreter','latex')
%% Upravljanje

t1 = zeros;
t2 = zeros;

Y1 = 1600;
v1u = 18;
v1i = 12.8;

Y2 = 1400;
v2u = 31;
v2i = 21;

P1 = zeros;
P2 = zeros;

dd1 = 0;
dd2 = 0;

% m2 = zeros;
% m1 = eye(55);

e1 = zeros;
e2 = zeros;

for k = 1:55
    
    if k == 55
        a = 1;
    else
        
    P1(k) = lc1(k+1) - lc1(k);
    P2(k) = lc2(k+1) - lc2(k);
    
    t1(k) = floor(P1(k)/v1i);
    t2(k) = floor(P2(k)/v2i);
    
    e1(k) = P1(k) - t1(k)*v1i;
    e2(k) = P2(k) - t2(k)*v2i;
    
        if dd1 > v1i
            t1(k) = t1(k) + floor(dd1/v1i);
            dd1 = dd1 - v1i*floor(dd1/v1i);
        else
            dd1 = e1(k) + dd1;
        end
        
        if dd2 > v2i
            t2(k) = t2(k) + floor(dd2/v2i);
            dd2 = dd2 - v2i*floor(dd2/v2i);
        else
            dd2 = e2(k) + dd2;
        end
    end
end
%% Pametni PWM

figure(8)
hold on

duty1 = zeros;
duty2 = zeros;
offset=0.5;
amp=0.5;
tm = zeros;
f1 = zeros;
f2 = zeros;
t0 = zeros;
tb = zeros;

for k = 1:length(t1)

    tm(k) = max(t1(k),t2(k));
    
    if k == length(t1)
        a = 2;    
    else
        if k == 1
            if t1(k) >= t2(k)
                
                f1(k) = 1/2;
                f2(k) = 1/2 * 1/(t1(k)-t2(k)+1);
                
                duty1(k) = 50;
                duty2(k) = 50 * 1/(t1(k)-t2(k)+1);
                
                if t2(k) == 0 
                    duty2(k) = 0;
                end
                
                t0(k) = 0;
                
                tb(k) = t0(k) + tm(k); 
                 
                t=t0(k)*2:0.0001:tb(k)*2;
        
                sq_wav1=offset+amp*square(2*pi*f1(k).*t,duty1(k));
                
                sq_wav2=offset+amp*square(2*pi*f2(k).*t,duty2(k));
                
                subplot(2,1,1)
                hold on
                plot(t,sq_wav1,'r')
                
                subplot(2,1,2)
                hold on
                plot(t,sq_wav2,'b')
                
                t0(k+1) = tb(k);
            end 
        else           
            if t1(k) > t2(k)
                
                f1(k) = 1/2;
                f2(k) = 1/2 * 1/(t1(k)-t2(k)+1);
                
                duty1(k) = 50;
                duty2(k) = 50 * 1/(t1(k)-t2(k)+1);
                
                if t2(k) == 0 
                    duty2(k) = 0;
                    f2(k) = f1(k);
                end
                
                tb(k) = t0(k) + tm(k); 
                 
                t=t0(k)*2:0.01:tb(k)*2;
        
                sq_wav1=offset+amp*square(2*pi*f1(k).*t,duty1(k));
                
                sq_wav2=offset+amp*square(2*pi*f2(k).*t,duty2(k));
                
                subplot(2,1,1)
                hold on
                plot(t,sq_wav1,'r')
                
                subplot(2,1,2)
                hold on
                plot(t,sq_wav2,'b')
                
                t0(k+1) = tb(k);
                
            else
                
                f2(k) = 1/2;
                f1(k) = 1/2 * 1/(t2(k)-t1(k)+1);
                
                duty2(k) = 50;
                duty1(k) = 50 * 1/(t2(k)-t1(k)+1);
                
                if t1(k) == 0 
                    duty1(k) = 0;
                end
                
                tb(k) = t0(k) + tm(k); 
                 
                t=t0(k)*2:0.01:tb(k)*2;
        
                sq_wav1=offset+amp*square(2*pi*f1(k).*t,duty1(k));
                
                sq_wav2=offset+amp*square(2*pi*f2(k).*t,duty2(k));
                
                subplot(2,1,1)
                axis([0 200 0 5])
                hold on
                plot(t,sq_wav1,'r')
                
                subplot(2,1,2)
                axis([0 200 0 5])
                hold on
                plot(t,sq_wav2,'b')
                
                t0(k+1) = tb(k);
            end
        end
    end
end
%% Ciscenje PWM Grafika

%Cilindar 1
subplot(2,1,1)
axis([0 200 0 5])
hold on
line([61.99 62],[0 1],'Color','white')      
line([191.99 192],[0 1],'Color','white')
line([65.99 66],[0 1],'Color','red')      

%Cilindar 2
subplot(2,1,2)
axis([0 200 0 5])
hold on

line([45.99 46],[0 1],'Color','white')
line([57.99 58],[0 1],'Color','white') 
line([58.99 59],[0 1],'Color','white')      
line([73.99 74],[0 1],'Color','white')      
line([85.99 86],[0 1],'Color','white')
line([95.99 96],[0 1],'Color','white')      
line([99.99 100],[0 1],'Color','white')      
line([105.99 106],[0 1],'Color','white')      
line([111.99 112],[0 1],'Color','white')      
line([131.99 132],[0 1],'Color','white')      
line([139.99 140],[0 1],'Color','white')
line([167.99 168],[0 1],'Color','white')      

line([49.99 50],[0 1],'Color','blue')
line([61.99 62],[0 1],'Color','blue')
line([77.99 78],[0 1],'Color','blue')
line([89.99 90],[0 1],'Color','blue')
line([97.99 98],[0 1],'Color','blue')
line([103.99 104],[0 1],'Color','blue')
line([109.99 110],[0 1],'Color','blue')
line([127.99 128],[0 1],'Color','blue')
line([135.99 136],[0 1],'Color','blue')
line([161.99 162],[0 1],'Color','blue')
%% Opterecenje geometrije masa na cilindru 1

g = 9.81;   %m/s2

mh = 300;    %kg, masa grabulje

m1 = 1602;   %kg
m2 = 1915.2; %kg, masa segmenta 2 sa grabuljom

mg = 500;   %kg, masa smeca

mc1 = 307.7; %kg, masa cilindra 1
mc2 = 160.2; %kg, masa cilindra 2

G1 = m1*g;       %N   Sila teze tela 1
G2 = m2*g;       %N   Sila teze tela 2
G2G = (m2+mg)*g; %N   Sila teze tela 2 sa grabuljom i smecem 
GG = mg*g;       %N   Sila teze smeca 

rho = 900; %kgm3

dk1 = 180;  %mm
ak1 = dk1^2/4*pi*0.000001;  %m2

dk2 = 140;  %mm
ak2 = dk2^2/4*pi*0.000001;  %m2

for k = 1:length(xH)

    alfa1(k) = pi/2 - atan(ygc1(k)/xgc1(k)); %cilin1
    alfa2(k) = pi/2 - atan(ygc2(k)/xgc2(k)); %cilin2
    alfa3(k) = pi/2 - atan(yg1(k)/xg1(k)); %c1
    alfa4(k) = pi/2 - atan(yg2(k)/xg2(k)); %c2
    alfa5(k) = atan((ygc1(k)-ya1(k))/(xa1(k)-xgc1(k))) + atan(ygc1(k)/xgc1(k)); %napadni ugao Fk, cilin1
    alfa6(k) = pi/2 - atan(yH(k)/xH(k));
    
    rg1(k,1) = xg1(k);
    rg1(k,2) = yg1(k);
    
    rg2(k,1) = xg2(k);
    rg2(k,2) = yg2(k);

    rgc1(k,1) = xgc1(k);
    rgc1(k,2) = ygc1(k);
        
    rgc2(k,1) = xgc2(k);
    rgc2(k,2) = ygc2(k);
    
    lg11(k) = sqrt(rg1(k,1).^2 + rg1(k,2).^2)*0.001;    %m
    lg22(k) = sqrt(rg2(k,1).^2 + rg2(k,2).^2)*0.001;    %m
    lg33(k) = sqrt(xH(k).^2 + yH(k).^2)*0.001; %m
    lgc11(k) = sqrt(rgc1(k,1).^2 + rgc1(k,2).^2)*0.001; %m
    lgc22(k) = sqrt(rgc2(k,1).^2 + rgc2(k,2).^2)*0.001; %m
    
    y11(k) = (lc1(k) - lc1(1))*0.001;  %m; trenutni hod cilindra 1
    y22(k) = (lc2(k) - lc2(1))*0.001;  %m; trenutni hod cilindra 2
    
    GC1(k) = (y11(k)*ak1*rho + mc1)*g;  %N
    GC2(k) = (y22(k)*ak2*rho + mc2)*g;  %N
        
    MGC1(k) = GC1(k)*lgc11(k)*sin(alfa1(k));    %Moment sile teze cilindra 1
    MGC2(k) = GC2(k)*lgc22(k)*sin(alfa2(k));    %Moment sile teze cilindra 2
    MG1(k) = G1*lg11(k)*sin(alfa3(k));          %Moment sile teze tela 1
    MG2(k) = G2*lg22(k)*sin(alfa4(k));          %Moment sile teze tela 2, + grabulja
    
    MG2G(k) = G2G*lg22(k)*sin(alfa4(k));        %Moment sile teze tela 2, smece + grabulja
    MG22(k) = (G2-mh*g)*lg22(k)*sin(alfa4(k));  %Moment sile teze tela 2
    MGG(k) = (GG+mh*g)*lg33(k)*sin(alfa6(k));   %Moment sile teze, grabulja + smece
    
    FLc1(k) = (MGC1(k) + MGC2(k) + MG1(k) + MG2G(k))/(lgc11(k)*sin(alfa5(k)));          %Sila na klipu,telo2 + smece + gr
    FLc2(k) = (MGC1(k) + MGC2(k) + MG1(k) + MG22(k) + MGG(k))/(lgc11(k)*sin(alfa5(k))); %Sila na klipu, poseban krak gr
    
    MLc1(k) = FLc1(k)*lgc11(k)*sin(alfa5(k));   %Moment sile, klip, telo2 + smece + grabulja
    MLc2(k) = FLc2(k)*lgc11(k)*sin(alfa5(k));   %Moment isle, klip, poseban krak smece + grabulja
    
end
    figure(15)
    
    subplot(3,1,1)
    hold on
    grid on
    plot(rad2deg(abs(q1-q1(1))),FLc1)
    title('Zavisnost $$F_L$$ cilindra 1 od $$\Delta \varphi_1$$ (uracunata masa tela 2)','interpreter','latex')
    xlabel('$$\Delta \varphi_1 [^{\circ}]$$','interpreter','latex')
    ylabel('$$F_L [N]$$','interpreter','latex')

    subplot(3,1,2)
    hold on
    grid on
    plot(rad2deg(abs(q1-q1(1))),MLc1)
    title('Zavisnost $$M_L$$ cilindra 1 od $$\Delta \varphi_1$$ (uracunata masa tela 2)','interpreter','latex')
    xlabel('$$\Delta \varphi_1 [^{\circ}]$$','interpreter','latex')
    ylabel('$$M_L [N]$$','interpreter','latex')
   
    subplot(3,1,3)
    hold on
    grid on
    plot(rad2deg(abs(q1-q1(1))),FLc1/ak1/100000)
    title('Zavisnost $$p_k$$ cilindra 1 od $$\Delta \varphi_1$$ (uracunata masa tela 2)','interpreter','latex')
    xlabel('$$\Delta \varphi_1 [^{\circ}]$$','interpreter','latex')
    ylabel('$$p_k [bar]$$','interpreter','latex')
    
    figure(16)
    subplot(3,1,1)
    hold on
    grid on
    plot(rad2deg(abs(q1-q1(1))),FLc2)
    title('Zavisnost $$F_L$$ cilindra 1 od $$\Delta \varphi_1$$ (sa potegom grabulje)','interpreter','latex')
    xlabel('$$\Delta \varphi_1 [^{\circ}]$$','interpreter','latex')
    ylabel('$$F_{L2} [N]$$','interpreter','latex')
    
    subplot(3,1,2)
    hold on
    grid on
    plot(rad2deg(abs(q1-q1(1))),MLc2)
    title('Zavisnost $$F_L$$ cilindra 1 od $$\Delta \varphi_1$$ (sa potegom grabulje)','interpreter','latex')
    xlabel('$$\Delta \varphi_1 [^{\circ}]$$','interpreter','latex')
    ylabel('$$M_{L2} [N]$$','interpreter','latex')
    
    subplot(3,1,3)
    hold on
    grid on
    plot(rad2deg(abs(q1-q1(1))),FLc2/ak1/100000)
    title('Zavisnost $$F_L$$ cilindra 1 od $$\Delta \varphi_1$$ (sa potegom grabulje)','interpreter','latex')
    xlabel('$$\Delta \varphi_1 [^{\circ}]$$','interpreter','latex')
    ylabel('$$p_{k2} [N]$$','interpreter','latex')
%% Sila reakcije resetke, trazenje polozaja
    
figure(20)
hold on
axis ([-15000 9000 -13500 5500])
% L = [];
for i = 50:55

    line([0, x1(50)], [0, y1(50)], 'LineWidth',2)
    line([0, x1(53)], [0, y1(53)], 'LineWidth',2)
    line([0, x1(55)], [0, y1(55)], 'LineWidth',2)
    line([x1(50), x2(55)], [y1(50), y2(55)], 'LineWidth',2)
    line([x1(53), x2(55)], [y1(53), y2(55)], 'LineWidth',2)
    line([x1(55), x2(55)], [y1(55), y2(55)], 'LineWidth',2)
    plot(0, 0, 'k^', 'MarkerSize',10, 'LineWidth',1.5)
    plot(x1(50), y1(50), 'kx', 'MarkerSize',10, 'LineWidth',1.5)
    plot(x1(53), y1(53), 'kx', 'MarkerSize',10, 'LineWidth',1.5)
    plot(x2(55), y2(55), 'rx', 'MarkerSize',10, 'LineWidth',1.5)
    plot(xH,yH,'k*')
    
    plot(xa2(50),ya2(50),'g*')
    plot(xa2(53),ya2(53),'g*')
    plot(xa2(55),ya2(55),'g*')
    
    plot(xa1(50),ya1(50),'g*')
    plot(xa1(53),ya1(53),'g*')
    plot(xa1(55),ya1(55),'g*')

    plot(xg1(50),yg1(50),'m*')
    plot(xg1(53),yg1(53),'m*')
    plot(xg1(55),yg1(55),'m*')
    
    plot(xg2(50),yg2(50),'m*')
    plot(xg2(53),yg2(53),'m*')
    plot(xg2(55),yg2(55),'m*')
    
    plot(xc2(50),yc2(50),'k*')
    plot(xc2(53),yc2(53),'k*')
    plot(xc2(55),yc2(55),'k*')
  
    kk(i) = (y1(i) - y2(55))/(x1(i) - x2(55));
    n(i) = -kk(i)*x2(55) + y2(55);
    phiii(i) = atan2((y1(i) - y2(55)),(x1(i) - x2(55)));

    xlin = -2891:0.1:-2650;
    yte = kk(i).*xlin + n(i);
    plot(xlin,yte)
    
%     for a = 50:55
%         for i = 1:1:length(xlin)
%     
%         L(i,a) = sqrt((x1(a) - xlin(i))^2+(y1(a) - yte(i))^2);
%     
%         if L(i,a) > 7077.5 && L(i,a) < 7077.7
%             AAA(i,a) = L(i,a);
%         end
%    
%         end
%     end

end

% AAAA = AAA(:,(50:55));


plot(C1(1), C1(2), 'k*', 'MarkerSize',8, 'LineWidth',1.5)

plot(xlin(2325), kk(50)*xlin(2325) + n(50), 'r^', 'MarkerSize',8, 'LineWidth',1.5)
% plot(xlin(1900), kk(51)*xlin(1900) + n(51), 'r^', 'MarkerSize',8, 'LineWidth',1.5)
% plot(xlin(1457), kk(52)*xlin(1457) + n(52), 'r^', 'MarkerSize',8, 'LineWidth',1.5)
plot(xlin(995), kk(53)*xlin(995) + n(53), 'r^', 'MarkerSize',8, 'LineWidth',1.5)
% plot(xlin(513), kk(54)*xlin(513) + n(54), 'r^', 'MarkerSize',8, 'LineWidth',1.5)
plot(xlin(1), kk(55)*xlin(1) + n(55), 'r^', 'MarkerSize',8, 'LineWidth',1.5)

line([xlin(1), xlin(2325)], [kk(55)*xlin(1) + n(55), kk(50)*xlin(2325) + n(50)], 'LineWidth',2)
% line([xlin(1), xlin(1900)], [kk(55)*xlin(1) + n(55), kk(51)*xlin(1900) + n(51)], 'LineWidth',2)
% line([xlin(1), xlin(1457)], [kk(55)*xlin(1) + n(55), kk(52)*xlin(1457) + n(52)], 'LineWidth',2)
line([xlin(1), xlin(995)], [kk(55)*xlin(1) + n(55), kk(53)*xlin(995) + n(53)], 'LineWidth',2)
% line([xlin(1), xlin(513)], [kk(55)*xlin(1) + n(55), kk(54)*xlin(513) + n(54)], 'LineWidth',2)


xlabel('$$x$$','interpreter','latex','fontsize',14)
ylabel('$$y$$','interpreter','latex','fontsize',14)
%% Prelivni pritisci za zaglavljeno stanje

FRn = 1500 * 2;  %150kg, sila reakcije resetke

for i = 50:55
    
    %telo 2 opterecenje
    
    rR(i) = sqrt((xH(55) - x1(i))^2 + (yH(55) - y1(i))^2);
    alfaL(i) = phiii(i) - pi/2; 
    FL22(i) = -(G2*(l2/2)*sin(alfaL(i)) + (GG+mh*g)*l2*sin(alfaL(i)) - FRn*rR(i))/(lb(i)*sin(pi - beta(i)));
    
    %sila na klipnjaci = moment teze tela 2 + moment smeca - sila reakcije resetke / krak sile na klipnjaci    
    
    pl22 = FL22./Ak;
    pL22 = pl22./1e5;

    %telo 1 opterecenje
    
    %alpha1
    %ugao izmedju za cilindar 1 i krak O1A1
    
    %y1(A1C1)
    kkc1a1(i) = (ya1(i) - C1(2))/(xa1(i) - C1(1));
    kkc1a1a(i) = atan2((ya1(i) - C1(2)),(xa1(i) - C1(1)));
    rad2deg(kkc1a1a);
    
    %y2(A1O1)
    kko1a1(i) = ya1(i)/xa1(i);
    kko1a1a(i) = atan2(ya1(i),xa1(i));
    rad2deg(kko1a1a);
   
    phic(i) = kkc1a1(i) - kko1a1(i);
    phici(i) = pi - atan2((kko1a1(i) - kkc1a1(i)),(1 + kko1a1(i)*kkc1a1(i)))
    rad2deg(phici);
    
    %alpha2
    %ugao za silu teze i krak O1MG1
    
    kkc1(i) = (yg1(i))/(xg1(i));
    kkc1a(i) = pi - (atan2(y1(i),x1(i)) - pi/2);
    rad2deg(kkc1a);
    
    %alpha3
    %ugao za krak O1A2
    
    kka2(i) = ya2(i)/xa2(i);
%     kkc2a(i) = atan2(ya2(i),xa2(i));
%     rad2deg(kkc2a)
    
    %ugao za silu clindra 2 i kraka O2A2
    kkc2(i) = (ya2(i) - yc2(i))/(xa2(i) - xc2(i));
%     kkc2c(i) = atan2((ya2(i) - yc2(i)),(xa2(i) - xc2(i)));
%     rad2deg(kkc2c)
    
    kaici(i) = pi - atan2((kkc2(i) - kka2(i)),(1 + kka2(i)*kkc2(i)));
    rad2deg(kaici)

    LBC(i) = sqrt(xa2(i).^2 + ya2(i).^2);
    LBC1(i) = sqrt(xa1(i).^2 + ya1(i).^2);

    FL11(i) = ( FL22(i)*LBC(i)*sin(kaici(i)) - G1*l1/2*sin(kkc1a(i)) ) / ( LBC1(i)*sin(phici(i)) ); 
    pl11 = FL11./ak1;
    pL11 = pl11./1e5;
end

figure(21)
hold on
grid on
title('Pritisak na strani klipnjache cilindra 2 (zaglavljen)','interpreter','latex')
plot(rad2deg(q2(50:end)-q2(1)),pL22(50:end))
xlabel('$$\varphi_2 [^{\circ}]$$','interpreter','latex')
ylabel('$$p_kp [bar]$$','interpreter','latex')

figure(22)
hold on
grid on
title('Pritisak na strani klipnjache cilindra 1 (zaglavljen)','interpreter','latex')
plot(rad2deg(q2(50:end)-q2(1)),pL11(50:end))
xlabel('$$\varphi_2 [^{\circ}]$$','interpreter','latex')
ylabel('$$p_kp [bar]$$','interpreter','latex')