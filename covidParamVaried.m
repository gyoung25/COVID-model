%Solves the ODE model found in Young et al 2021 (F1000Research) with parameters varied up to 10% of baseline, simulates the model 1000 times and plots the baseline infected-class solution curve along with two curves that are representative of 80% of all simulations

function covidParamVaried()

%params
phiVec=linspace(0,1,6);
xi=.94;
mu=0.02;
a=0;
c=0.85;
%params=[beta,sigma,rho,gammaA,gammaM,gammaH,gammaC,delta1,delta2,delta3,m,x1,x2,x3];
%baseline values
BLParams = [0.76, 0.32, 0.37, .23, .23, .22, .23, 0.69, .15, .55,...
    .17, .18, .25, .51];

N=1000;%number of simulations with random params
IcMax=zeros(N,1);
tRef=linspace(0,100,5000);
IcInterpHigh=zeros(N,5000);
IcInterpLow=zeros(N,5000);
% figure;
% hold on

for j=1:length(phiVec)
        phi=phiVec(j);
    for i=1:N
        %on the first iteration, use baseline params
        if i==1
            paramVals=BLParams;
            params = num2cell(paramVals);
            [beta, sigma, rho, gammaA, gammaM, gammaH, gammaC, delta1, delta2, delta3,...
            m, x1, x2, x3]=deal(params{:});
        else
            %vary each param value randomly by up to +/-(rRang)%
            pRang=0.1;
            paramVals=BLParams+(-pRang+2*pRang*rand(size(BLParams))).*BLParams;
            params = num2cell(paramVals);
            [beta, sigma, rho, gammaA, gammaM, gammaH, gammaC, delta1, delta2, delta3,...
            m, x1, x2, x3]=deal(params{:});
        end

        %the time Matlab will solve the ODE until
        tMax=100;

        %initial value of y
        %[S X V E Ia Ip Im Ih Ic D2 R2]
        y0=[0.99 0 0 0.01 0 0 0 0 0 0 0];

        options = odeset('RelTol',1e-10,'AbsTol',1e-12);
        [t,y]=ode45(@f,[0 tMax],y0,options);

        Ic=y(:,9);

        IcMax(i)=max(Ic);
        IcInterp(i,:)=interp1(t,Ic,tRef);
        if i==1
            %plot(t,Ic,'b','linewidth',4)
            %length(Ic)
            IcRef=interp1(t,Ic,tRef);
        else
            %plot(t,Ic,'linewidth',2)
            if IcMax(i)>IcMax(1)
                IcInterpHigh(i,:)=IcInterpHigh(i,:)+IcInterp(i,:);
            else
                IcInterpLow(i,:)=IcInterpLow(i,:)+IcInterp(i,:);
            end

        end


    end
    [IcMax sortInd]=sort(IcMax);
    ICHigh(j,:)=IcInterp(sortInd(round(N*0.9)),:);
    ICLow(j,:)=IcInterp(sortInd(round(N*0.1)),:);
    ICOG(j,:)=IcInterp(sortInd(find(sortInd==1)),:);
end
%Plot the results: plot(x,y) produces a figure of array x versus array y,
%both of which must be the same length. You can also include options in the
%plot function, like color, linewidth (thickness), etc
% figure;
% plot(t,Ic,'linewidth',2)
% hold on %do not start new figure
%plot(t,yExact,'r*','linewidth',2)
% xlabel('time (days)','fontsize',18)
% ylabel('ICU cases','fontsize',18)
% set(gca,'fontsize',20) %this sets the axis label font size


% std(IcMax)/IcMax(1)
% oneStdDev=sum(IcMax<IcMax(1)+std(IcMax))/N
% within22=sum(IcMax<1.22*IcMax(1))/N
% within30=sum(IcMax<1.3*IcMax(1))/N
% within40=sum(IcMax<1.4*IcMax(1))/N

% aveHigh=sum(IcMax(find(IcMax>IcMax(1))))/sum(IcMax>IcMax(1))
% aveLow=sum(IcMax(find(IcMax<IcMax(1))))/sum(IcMax<IcMax(1))
% countHigh=sum(IcMax>IcMax(1))
% countLow=sum(IcMax<IcMax(1))
% IcMax(1)

IcInterpHigh=IcInterpHigh(any(IcInterpHigh,2),:);
IcInterpLow=IcInterpLow(any(IcInterpLow,2),:);

IcHighAve=mean(IcInterpHigh);
IcLowAve=mean(IcInterpLow);
% 
% figure
% hold on
% plot(tRef,IcHighAve,'linewidth',2)
% plot(tRef,IcRef,'linewidth',2)
% plot(tRef,IcLowAve,'linewidth',2)
% xlabel('time (days)','fontsize',18)
% ylabel('ICU cases','fontsize',18)
% set(gca,'fontsize',20)


orderedIcMax=sort(IcMax);
[IcHist edges]=histcounts(orderedIcMax,'BinWidth',0.0005)
IcPMF=IcHist/N;
IcCDF=cumsum(IcPMF)
%the next two lines find the indeces between which 80% of IcMaxes occur
[min10 iMinHist]=min(abs(IcCDF-0.1));
[min90 iMaxHist]=min(abs(IcCDF-0.9));
[dum iMinIC]=min(abs(orderedIcMax-edges(iMinHist)));
[dum iMaxIC]=min(abs(orderedIcMax-edges(iMaxHist)));
IcMax80=orderedIcMax(iMinIC+1:iMaxIC-1);
% figure
% histogram(IcMax,'BinWidth',0.0005)
% hold on
% histogram(IcMax80,'BinWidth',0.0005)
%plot([edges(iMinHist) edges(iMinHist)],[0 IcHist(iMinHist)],'r','linewidth',4)

A4=load('timeSeries_phiPt6a0.dat');
t4=A4(:,1);
E4=A4(:,5);
Ic4=A4(:,10);

fig=figure
hold on
plot(t4,Ic4,'k','linewidth',2)
%plot(tRef,ICOG(4,:),'b--','linewidth',2)
plot(tRef,ICHigh(4,:),'k--','linewidth',2)
plot(tRef,ICLow(4,:),'k--','linewidth',2)
axis([0 100 0 0.008])
ylabel('ICU cases','fontsize',18)
xlabel('time (days)','fontsize',18)
set(gca,'fontsize',20)
%legend('\phi=0, a=0','\phi=0.6, a=0','\phi=0, a=0.8','\phi=0.6, a=0.8','location','northeast','fontsize',23)
box on
%print(fig,'-depsc','sim80prcRange.eps')

%A0=load('ICU_peak_phiVSIC_a0_xiPt94_cPt85.dat');
IcMaxLow=max(ICLow,[],2);
IcMaxHigh=max(ICHigh,[],2);
IcMaxOG=max(ICOG,[],2)
fig=figure;
hold on
%plot(A0(:,1),A0(:,2),'r','linewidth',5)
plot(phiVec,IcMaxOG,'k','linewidth',2)
plot(phiVec,IcMaxLow,'o-k','markerfacecolor','k','linewidth',2)
plot(phiVec,IcMaxHigh,'o-k','markerfacecolor','k','linewidth',2)
axis([0 1 0 0.01])
ylabel('Peak ICU','fontsize',18)
xlabel('\phi','fontsize',18)
set(gca,'fontsize',20)
set(gca,'YScale','log')
%legend('a=0', 'a=0.5','a=0.8','location','southeast','fontsize',23)
box on
%print(fig,'-depsc','ICUpeaks80rcRange_logScale.eps')

    % yPrime is the variable returned by the function
    % f is the name of the function
    % t and y are the inputs
    function yPrime=f(t,y)
        S=y(1);
        X=y(2);
        V=y(3);
        E=y(4);
        Ia=y(5);
        Ip=y(6);
        Im=y(7);
        Ih=y(8);
        Ic=y(9);
        D1=y(10);
        R2=y(11);
        
        Sp=-(1-a*c)*beta*S*(Ia+Ip+Im)-mu*S;
        Xp=mu*(1-phi)*S-(1-a*c)*beta*X*(Ia+Ip+Im);
        Vp=mu*phi*S-(1-xi)*(1-a*c)*beta*V*(Ia+Ip+Im);
        Ep=(1-a*c)*beta*(Ia+Ip+Im)*(S+X+(1-xi)*V)-sigma*E;
        Iap=sigma*rho*E-gammaA*Ia;
        Ipp=sigma*(1-rho)*E-delta1*Ip;
        Imp=delta1*Ip-x1*delta2*Im-(1-x1)*gammaM*Im;
        Ihp=x1*delta2*Im-x2*delta3*Ih-(1-x2)*gammaH*Ih;
        Icp=x2*delta3*Ih-(1-x3)*gammaC*Ic-x3*m*Ic;
        D1p=x3*m*Ic;
        R2p=gammaA*Ia+(1-x1)*gammaM*Im+(1-x2)*gammaH*Ih+(1-x3)*gammaC*Ic;
        
        yPrime=[Sp Xp Vp Ep Iap Ipp Imp Ihp Icp D1p R2p]';
    end


end