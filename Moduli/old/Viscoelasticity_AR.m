tic
% clear_all_but NMSD_x NMSD_y A pName fNames;
delete output1.mat;
delete output2.mat;
delete output3.mat;
delete output4.mat;
%lambda=(1*10^-6)/(6*pi*2.5*10^-6*0.896*10^-3);
lambda=(1.38064852*(10^-23)*293.15)/(pi*0.5*10^-6);
%AR=33;%Aquisition rate
AT=MSD(end,1);%Recorded time
AR=numel(MSD(:,1))/AT;%Aquisition rate
tau=1/AR:1/AR:AT;
subplot(3,1,1);
loglog(tau,MSD(:,2)*10^-12,'go');
AA=MSD(:,2)*10^-12;
MSD_smooth=1;
AA=smoothdata(AA,'gaussian',MSD_smooth);
Ceil_FN=1.15;
jj=2;
bb2(1,1)=1/AR;
bb2(1,2)=AA(1);
bb=zeros(1,(numel(tau)));
for pp=2:numel(tau)

    bb(pp)=ceil(Ceil_FN^(pp));%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if bb(pp)~=bb(pp-1)
    if bb(pp)<=numel(tau)
        bb2(jj,1)=bb(pp)/AR;
        bb2(jj,2)=AA(bb(pp));
         jj=jj+1;
    end
    end
    
end
%omega=logspace(-4,5,100);
%loglog(bb1(:,1),bb1(:,2),'ko');
hold on    
%subplot(2,1,1);
%bb1(:,2)=smoothdata(bb1(:,2),'gaussian',50);
%plot(bb1(:,1),bb1(:,2),'ro');
% hold on
%bb2(:,2)=smoothdata(bb2(:,2),'gaussian',NMSD_smooth);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%loglog(bb2(:,1),bb2(:,2),'ko');
loglog(bb2(:,1),bb2(:,2),'bo');
hold on

fs=1*10^5;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xx(:)=bb2(1,1):1/fs:bb2(end,1);
cs(:)=pchip(bb2(:,1),bb2(:,2),xx);
loglog(xx,cs,'b-');
%plot(xx,cs,'b-');
ylabel('\pi(\tau)','fontweight','bold','fontsize',16)
xlabel('\tau (sec)','fontweight','bold','fontsize',16)
hold on
jj=2;
bb3(1,1)=1/AR;
bb3(1,2)=cs(1);
clear bb;
bb=zeros(1,(numel(xx)));
for pp=2:numel(xx)

    bb(pp)=ceil(1.015^(pp));
    if bb(pp)~=bb(pp-1)
    if bb(pp)<=numel(xx)
        bb3(jj,1)=(1/AR)+((bb(pp)-1)/fs);
        %bb3(jj,1)=bb(pp)/fs;
        bb3(jj,2)=cs(bb(pp));
         jj=jj+1;
    end
    end
    
end
%loglog(bb3(:,1),bb3(:,2),'ro');
%ginf=((bb2(end,2)-bb2(end-9,2))/(bb2(end,1)-bb2(end-9,1)));
ginf=((bb3(end,2)-bb3(end-20,2))/(bb3(end,1)-bb3(end-20,1)));
%ginf=((bb3(end,2)-bb3(end-9,2))/(bb3(end,1)-bb3(end-9,1)));
%ginf=0;
omega=logspace(-4,5,500);
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parfor mv =1:4
    if mv==1
        cs1=cs;
        xx1=xx;
        g1=zeros(1,125);
        omega1=omega(1:125);        
        g1inf=ginf;
            
        for ss1=1:numel(omega1)
                disp(num2str(ss1));
            I1=zeros(1,(numel(xx1)-1));
            for k1=2:numel(xx1)
                I1(k1-1)=((cs1(k1)-cs1(k1-1))/(xx1(k1)-xx1(k1-1)))*(exp(-1i*omega1(ss1)*xx1(k1-1))-exp(-1i*omega1(ss1)*xx1(k1)));
            end
            g1(ss1)=(1i*omega1(ss1)*0)+((1-(exp(-1i*omega1(ss1)*xx1(1))))*((cs1(1)-0)/xx1(1)))+(g1inf*exp(-1i*omega1(ss1)*xx1(end)))+sum(I1);
            g1(ss1)=-g1(ss1)/(omega1(ss1)^2);
        end
        m=matfile(sprintf('output%d.mat', mv),'writable',true);
        m.g_1=g1;
        m.omega1=omega1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if mv==2
        cs2=cs;
        xx2=xx;
        g2=zeros(1,125);
        omega2=omega(126:250);        
        g2inf=ginf;
        for ss2=1:numel(omega2)
            I2=zeros(1,(numel(xx2)-1));
            for k2=2:numel(xx2)
                I2(k2-1)=((cs2(k2)-cs2(k2-1))/(xx2(k2)-xx2(k2-1)))*(exp(-1i*omega2(ss2)*xx2(k2-1))-exp(-1i*omega2(ss2)*xx2(k2)));
            end
            g2(ss2)=(1i*omega2(ss2)*0)+((1-(exp(-1i*omega2(ss2)*xx2(1))))*((cs2(1)-0)/xx2(1)))+(g2inf*exp(-1i*omega2(ss2)*xx2(end)))+sum(I2);
            g2(ss2)=-g2(ss2)/(omega2(ss2)^2);
        end
        m=matfile(sprintf('output%d.mat', mv),'writable',true);
        m.g_2=g2;
        m.omega2=omega2;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if mv==3
        cs3=cs;
        xx3=xx;
        g3=zeros(1,125);
        omega3=omega(251:375);        
        g3inf=ginf;
        for ss3=1:numel(omega3)
            I3=zeros(1,(numel(xx3)-1));
            for k3=2:numel(xx3)
                I3(k3-1)=((cs3(k3)-cs3(k3-1))/(xx3(k3)-xx3(k3-1)))*(exp(-1i*omega3(ss3)*xx3(k3-1))-exp(-1i*omega3(ss3)*xx3(k3)));
            end
            g3(ss3)=(1i*omega3(ss3)*0)+((1-(exp(-1i*omega3(ss3)*xx3(1))))*((cs3(1)-0)/xx3(1)))+(g3inf*exp(-1i*omega3(ss3)*xx3(end)))+sum(I3);
            g3(ss3)=-g3(ss3)/(omega3(ss3)^2);
        end
        m=matfile(sprintf('output%d.mat', mv),'writable',true);
        m.g_3=g3;
        m.omega3=omega3;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if mv==4
        cs4=cs;
        xx4=xx;
        g4=zeros(1,125);
        omega4=omega(376:500);        
        g4inf=ginf;
        for ss4=1:numel(omega4)
            I4=zeros(1,(numel(xx4)-1));
            for k4=2:numel(xx4)
                I4(k4-1)=((cs4(k4)-cs4(k4-1))/(xx4(k4)-xx4(k4-1)))*(exp(-1i*omega4(ss4)*xx4(k4-1))-exp(-1i*omega4(ss4)*xx4(k4)));
            end
            g4(ss4)=(1i*omega4(ss4)*0)+((1-(exp(-1i*omega4(ss4)*xx4(1))))*((cs4(1)-0)/xx4(1)))+(g4inf*exp(-1i*omega4(ss4)*xx4(end)))+sum(I4);
            g4(ss4)=-g4(ss4)/(omega4(ss4)^2);
        end
        m=matfile(sprintf('output%d.mat', mv),'writable',true);
        m.g_4=g4;
        m.omega4=omega4;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load output1.mat;
load output2.mat;
load output3.mat;
load output4.mat;

g_final=horzcat(g_1,g_2,g_3,g_4);
for q=1:numel(omega)
    gg(q)=((1/(1i*omega(q)*g_final(q))));
end

% subplot(4,1,2);

%loglog(omega,real(gg),'ro-');
hold on;
% loglog(omega,imag(gg),'b*-');
new1=real(gg);
new2=imag(gg);
%subplot(4,1,3);
Elastic_modulus=(new1*lambda)';
Viscous_modulus=(new2*lambda)';
% real_smooth=smoothdata(Elastic_modulus(202:300),'sgolay',10);
% imag_smooth=smoothdata(Viscous_modulus(202:300),'sgolay',10);
for i=1:numel(omega)
Complex_viscosity(i)=(sqrt(Elastic_modulus(i)^2+Viscous_modulus(i)^2))/omega(i);
end
%Zero_shear_viscosity=mean(Complex_viscosity(1:56));
%Zero_shear_viscosity_stderror=std(Complex_viscosity(1:56))/sqrt(length(Complex_viscosity(1:56)));
omega=omega';
subplot(3,1,2);
loglog(omega,Elastic_modulus,'ro-');
hold on;
loglog(omega,Viscous_modulus,'b*-');
subplot(3,1,3);
loglog(omega,Complex_viscosity,'ro-');
toc