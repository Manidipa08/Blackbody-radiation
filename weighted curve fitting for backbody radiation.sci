Temp = [2200 3000 4500 5000 5500 5580 6000 6500 7000 7800 8500]
xx=(1./Temp(1:5)); yy=w(cors_w(1:5)); W=ones(1,5)
XY=xx.*yy
XX=xx.^2
m=length(xx)
s_x=sum(xx)
s_y=sum(yy)
s_xy=sum(XY)
s_xx=sum(XX)
s_w=sum(W)
s_wx=sum(W.*xx)
s_wy=sum(W.*yy)
s_wxx=sum(W.*xx.*xx)
s_wxy=sum(W.*xx.*yy)
a0=((s_wy.*s_wxx)-(s_wxy.*s_wx))/((s_w.*s_wxx)-(s_wx)^2)//a0-intercept
a1=((s_w.*s_wxy)-(s_wx.*s_wy))/((s_w.*s_wxx)-(s_wx)^2)///a1-slope
[a b]=reglin(xx,yy)
fx=a0+a1*xx
