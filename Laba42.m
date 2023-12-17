clear all;
close all;
X=load('R4_08.txt');
S=length(X);
s1=X(1:S,1);
s2=X(1:S,2);
s3=X(1:S,3);

Fd=250;
tmax=S/Fd;
T=1/Fd;
t=0:T:tmax-T;
x0=20; 
y0=50;
w0=1800;
h0=900;
figure('position',[x0,y0,w0,h0])
x1=30;
dx=50;
dy=50;
y1=30;
w=650;
h=250;
hAxes11=axes('units','pixels','position',[x1,y1,w/2,h]);
hAxes21=axes('units','pixels','position',[x1+dx+w/2,y1,w,h]);
hAxes31=axes('units','pixels','position',[x1+2*dx+3/2*w,y1,w,h]);
hAxes12=axes('units','pixels','position',[x1,y1+dy+h,w/2,h]);
hAxes22=axes('units','pixels','position',[x1+dx+w/2,y1+dy+h,w,h]);
hAxes32=axes('units','pixels','position',[x1+2*dx+3/2*w,y1+dy+h,w,h]);
hAxes13=axes('units','pixels','position',[x1,y1+2*dy+2*h,w/2,h]);
hAxes23=axes('units','pixels','position',[x1+dx+w/2,y1+2*dy+2*h,w,h]);
hAxes33=axes('units','pixels','position',[x1+2*dx+3/2*w,y1+2*dy+2*h,w,h]);
axes(hAxes13)
plot(t,s1);
grid on;
hold on
m1=mean(s1);
s01=s1-m1;
w=hamming(S); % —оздание массива значений оконной функции УwФ 
for j=1:S 	% ”множение сигнала на оконную функцию: 
    sw1(j)=s01(j)*w(j); 
end 
plot(t,sw1)
ft1=fft(s1);
N=tmax*Fd;
for j=1:N
    if (j==1)
        as1(j)=sqrt(real(ft1(j))^2+imag(ft1(j))^2)/N;
    else
        as1(j)=sqrt(real(ft1(j))^2+imag(ft1(j))^2)/N*2;
    end
end
axes(hAxes23)
df=Fd/N;
for j=1:(N/2)
    f(j)=df*(j-1);
end
stem(f,as1(1:(N/2)),'.')
xlim([0,20]);
for j=1:N
    if (j==1)
        psd1(j)=1/N*(real(ft1(j))^2+imag(ft1(j))^2)/Fd;
    else
        psd1(j)=2/N*(real(ft1(j))^2+imag(ft1(j))^2)/Fd;
    end
end
axes(hAxes33)
plot(f,psd1(1:(N/2)))
xlim([0,20]);

axes(hAxes12)
plot(t,s2);
grid on
hold on
m2=mean(s2);
s02=s2-m2;
for j=1:S 	% ”множение сигнала на оконную функцию: 
    sw2(j)=s02(j)*w(j); 
end 
plot(t,sw2)
ft2=fft(s2);
for j=1:N
    if (j==1)
        as2(j)=sqrt(real(ft2(j))^2+imag(ft2(j))^2)/N;
    else
        as2(j)=sqrt(real(ft2(j))^2+imag(ft2(j))^2)/N*2;
    end
end
axes(hAxes22)
stem(f,as2(1:(N/2)),'.')
xlim([0,20]);
for j=1:N
    if (j==1)
        psd2(j)=1/N*(real(ft2(j))^2+imag(ft2(j))^2)/Fd;
    else
        psd2(j)=2/N*(real(ft2(j))^2+imag(ft2(j))^2)/Fd;
    end
end
axes(hAxes32)
plot(f,psd2(1:(N/2)))
xlim([0,20]);

axes(hAxes11)
plot(t,s3)
grid on
hold on
m3=mean(s3);
s03=s3-m3;
for j=1:S 	% ”множение сигнала на оконную функцию: 
    sw3(j)=s03(j)*w(j); 
end 
plot(t,sw3)
ft3=fft(s3);
for j=1:N
    if (j==1)
        as3(j)=sqrt(real(ft3(j))^2+imag(ft3(j))^2)/N;
    else
        as3(j)=sqrt(real(ft3(j))^2+imag(ft3(j))^2)/N*2;
    end
end
axes(hAxes21)
stem(f,as3(1:(N/2)),'.')
xlim([0,20]);
for j=1:N
    if (j==1)
        psd3(j)=1/N*(real(ft3(j))^2+imag(ft3(j))^2)/Fd;
    else
        psd3(j)=2/N*(real(ft3(j))^2+imag(ft3(j))^2)/Fd;
    end
end
axes(hAxes31)
plot(f,psd3(1:(N/2)))
xlim([0,20]);
periodogram
