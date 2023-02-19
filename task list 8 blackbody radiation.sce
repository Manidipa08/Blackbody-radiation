clc
clear
n=2000
v=linspace(0,1.2*10^14,n)
h=6.63D-34//plank constant
c=3*10^8
k=1.38*10^(-23)
T=500
//frequency combined
for i = 1:n
    E_P(i)=(8*%pi*h*(v(i)^3))/(c^3*(exp((h*v(i))/(k*T))-1))
    E_RJ(i)=(8*%pi*v(i)^2*k*T)/(c^3)
    E_W(i)=(8*%pi*h*(v(i)^3)*exp(-(h*v(i))/(k*T)))/(c^3)
end
show_window(1)
plot(v,E_P,'k','linewidth',2)
plot(v,E_RJ,'-*','linewidth',2)
plot(v,E_W,'.-r','linewidth',2)
m=gca();
m.data_bounds=[0 0;1.48*10^14 1.48*10^(-18)]
title("Energy Plots against frequency")
xgrid(7)
//wavelength combined
w=linspace(0,50*10^(-6),n)
for i=1:n
    E_PW(i)=(8*%pi*h*c)/(w(i)^5*(exp((h*c)/(w(i)*k*T))-1))
    E_RJW(i)=(8*%pi*k*T)/w(i)^4
    E_WW(i)=(8*%pi*h*c*exp(-(h*c)/(w(i)*k*T)))/w(i)^5
end
show_window(2)
plot(w,E_PW,'k','linewidth',2)
plot(w,E_RJW,'-*','linewidth',2)
plot(w,E_WW,'.-r','linewidth',1)
at=gca();
at.data_bounds=[2*10^(-6) 0;50*10^(-6) 1*7]
title("Energy Plots against wavelength")
xgrid(7)
//Plank for different temperature comparison against wavelength--------------------------
Temp = [2200 3000 4500 5000 5500 5580 6000 6500 7000 7800 8500]
for p = 1:11
    for i=1:n
        E_Planck(i,p)=(8*%pi*h*c)/(w(i)^5*(exp((h*c)/(w(i)*k*Temp(p)))-1))
        E_RelJ(i,p)=(8*%pi*k*Temp(p))/w(i)^4
        E_Wein(i,p)=(8*%pi*h*c*exp(-(h*c)/(w(i)*k*Temp(p))))/w(i)^5
    end
end
show_window(3)
plot(w,E_Planck(:,1:6),'-*','thickness',1.5); 
at=gca(); 
at.data_bounds=[0,0;3e-6,1e06]
xgrid(7)
for g=1:6
    legend([string(Temp(g))])
end
title("Energy_Planck Plots against wavelength at diff temp")
show_window(4)
plot(w,E_RelJ(:,1:6),'-*','thickness',1.5); 
at=gca(); 
at.data_bounds=[0,0;3e-6,8e06]
xgrid(7)
for g=1:6
    legend([string(Temp(g))])
end
title("Energy_RelJ Plots against wavelength at diff temp")
show_window(5)
plot(w,E_Wein(:,1:6),'-*','thickness',1.5); 
at=gca(); 
at.data_bounds=[0,0;3e-6,1e06]
xgrid(7)
for g=1:6
    legend([string(Temp(g))])
end
title("Energy_Wein Plots against wavelength at diff temp")
//Plank for different temperature comparison against frequency--------------------------
f=linspace(0,20*10^14,n)
for p = 1:11
    for i=1:n
        E_Planckf(i,p)=(8*%pi*h*(f(i)^3))/(c^3*(exp((h*f(i))/(k*Temp(p)))-1))
        E_RelJf(i,p)=(8*%pi*f(i)^2*k*Temp(p))/(c^3)
        E_Weinf(i,p)=(8*%pi*h*(f(i)^3)*exp(-(h*f(i))/(k*Temp(p))))/(c^3)
    end
end
show_window(6)
plot(f,E_Planckf(:,1:6),'-*','thickness',1.5); 
at=gca(); 
at.data_bounds=[0 0;2*10^15 1.4*10^(-15)]
title("Energy_Planck Plots against frequency at diff temp")
xgrid(7)
for g=1:6
    //h(g)=gce()
    legend([string(Temp(g))])
end
show_window(7)
plot(f,E_RelJf(:,1:6),'-*','thickness',1.5); 
at=gca(); 
at.data_bounds=[0 0;1.8*10^13 2.2*10^(-17)]
title("Energy_RelJ Plots against frequency at diff temp")
xgrid(7)
for g=1:6
    legend([string(Temp(g))])
end
show_window(8)
plot(f,E_Weinf(:,1:6),'thickness',1.5) 
at=gca(); 
at.data_bounds=[0 0;1.7*10^15 1.4*10^(-15)]
title("Energy_Wein Plots against frequency at diff temp")
xgrid(7)
for g=1:6
    legend([string(Temp(g))])
end
//Wien's Displacement law & Wien's Constant ---------------------
for y=1:11
    [x(y),cors_w(y)]=max(E_Planck(:,y))
    Wein_Cons(y)=w(cors_w(y))*Temp(y)
end
disp("Wein constant calculated for "+string(Temp'(1:4))+" K =  "+string(Wein_Cons(1:4)))
//least square fitting -- WEIN CONSTANT 
show_window(9)
xx=(1./Temp(1:5)); yy=w'(cors_w(1:5)); W=ones(1,5)

exec('C:\Users\MANIDIPA BANERJEE\Desktop\MP III Scilab SEM - IV\weighted curve fitting for backbody radiation.sci', -1)//exec('C:\Users\MANIDIPA BANERJEE\Desktop\MP III Scilab SEM - IV\weighted curve fitting for backbody radiation.sci', -1)
plot(xx,yy,'or')
plot(xx,fx,'*-k')
title('<<<< Weighted Least Square Fitting >>>>','color','Blue','Fontsize','5')
xlabel("X-values---------->",'color','green','Fontsize','4')
ylabel("Y-values---------->",'color','green','Fontsize','4')
legend('calculated','reglin')
disp("Wein Constant (using weighted least square fitting) = "+string(a1))
disp("Wein Constant (using reglin) = "+string(a))
//least square fitting --STEPHAN CONSTANT 
show_window(10)
xxx=(Temp(1:10).^4); 
for q=1:10
    yyy(q)=(c/4)*inttrap(w(2:n),E_Planck(:,q)(2:n)); 
end
W=ones(1,10)
exec('C:\Users\MANIDIPA BANERJEE\Desktop\MP III LAB IV files\Weighted lsf for stephan.sci', -1)
plot(xxx,yyy,'o-r','linewidth',2)
plot(xxx,fx,'*-y')
title('<<<< Weighted Least Square Fitting >>>>','color','Blue','Fontsize','5')
xlabel("X-values---------->",'color','green','Fontsize','4')
ylabel("Y-values---------->",'color','green','Fontsize','4')
legend('calculated','reglin')
disp("Stephan Constant (using weighted least square fitting) = "+string(a1))
disp("Stephan Constant (using reglin) = "+string(a))
