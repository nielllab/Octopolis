%%% code to analyze cephalopod optomotor response
%%% needs DLC tracked points and BORIS timing

clear all

d = readtable('PXL_20231003_235458873DLC_resnet50_squid_okr_100223Oct5shuffle1_100000.csv');

%%% a few parameters
fps = 30;   %%% video frames per second
pthresh = 0.97;   %%% likelihood threshold for DLC points

%% check out head tracking first, to get a sense of the data
headxy = d{:,20:21};  %%% pull data from DLC
headp = d{:,22};
figure
plot(headxy(:,1),headxy(:,2))

figure
plot(headp)
ylabel('head likelihood')
xlabel('frame')

%%%  threshold points on likelihood
headxy(headp<pthresh,:) = NaN;
figure
%plot(headxy(:,1),headxy(:,2))
plot(headxy(:,1))

%%% interpolate over nans
headxy(:,1) = interpNan(headxy(:,1),20);
figure
plot(headxy(:,1))

%% analyse head angle based on axis between the two eyes

%%% load/filter left eye data
eyeL = d{:,2:3};
eyeLp = d{:,4};
eyeL(eyeLp<pthresh,:) = NaN;
for i = 1:2
    eyeL(:,i) = interpNan(eyeL(:,i),300,'pchip');
end
figure
plot(eyeL(:,1),eyeL(:,2));

%%% load/filter right eye data
eyeR = d{:,14:15};
eyeRp = d{:,16};
eyeR(eyeLp<pthresh,:) = NaN;
for i = 1:2
    eyeR(:,i) = interpNan(eyeR(:,i),300);
end
figure
plot(eyeR(:,1),eyeR(:,2)); hold on
plot(eyeL(:,1),eyeL(:,2));

%%% calculate angle between the two eyes
th = atan2(eyeR(:,2)-eyeL(:,2), eyeR(:,1) - eyeL(:,1));
th = th*180/pi;

t = (1:length(th))*(1/30);
figure
plot(t,th)
xlabel('secs'); ylabel('theta')
xlim([0 120])

%%% calculate speed based on change in angle over time interval (defined by 'gap')
%dth = diff(th);
gap =10;
dth = th((gap+1):end) - th(1:(end-gap));
dth(dth>90) = dth(dth>90)-360;
dth(dth<-90) = dth(dth<-90)+360;
dth = dth/(gap/30);
figure
plot(dth)

t = (1:length(dth))*(1/30);
dth_smooth = medfilt1(dth,31);
figure
plot(t,dth); hold on
plot(t,dth_smooth,'Linewidth',2)
axis([0 120 -100 100])

%% read in Boris events
b = readtable('8873VideoEvents.csv');
eventID = b{:,11};
rt = 10;
eventF = b{:,18}-rt; %%% boris labeling reaction time
eventF = eventF; 
cwstartF = eventF(eventID==1);
ccwstartF = eventF(eventID==2);
stopF = eventF(eventID==3);

%%% plot head speed with Boris timing
figure
t = (1:length(dth_smooth))/fps;
plot(t,dth_smooth,'Linewidth',2); hold on;
for i = 1:length(cwstartF)
    plot(cwstartF(i)/30,0,'b*')
end
for i = 1:length(ccwstartF)
    plot(ccwstartF(i)/30,0,'g*')
end
for i = 1:length(stopF)
    plot(stopF(i)/30,0,'r*')
end
axis([0 120 -100 100])


%%% use Boris events to define a speed variable
eventF = eventF(~isnan(eventID));
eventID = eventID(~isnan(eventID));
stim = zeros(size(dth));
for i = 1:length(eventID)-1;
    stim(eventF(i):eventF(i+1)) = eventID(i);
end
figure
plot(stim)
stimV=zeros(size(stim));
stimV(stim==1)=50;
stimV(stim==2) = -50;
stimV(stim==3) = 0;
figure
plot(stimV)

figure
plot(t,dth_smooth);
hold on
plot(t,stimV);

%%% calculate speed in different conditions
nframes = 300*fps   %%% option to analyze a subset of data (for grant, used first 5 mins which were one SF)
dth_short = dth_smooth(1:nframes);
stim_short = stimV(1:nframes);
spd(1) = nanmean(dth_short(stim_short==50));
spd(2) = nanmean(dth_short(stim_short==-50));
spd(3) = nanmean(dth_short(stim_short==0));

nstim = length(find(eventF<nframes & eventID==1))

spd_std(1) = nanstd(dth_short(stim_short==50))/sqrt(length(find(eventF<nframes & eventID==1)));
spd_std(2) = nanstd(dth_short(stim_short==-50))/sqrt(length(find(eventF<nframes & eventID==2)));
spd_std(3) = nanstd(dth_short(stim_short==0))/sqrt(length(find(eventF<nframes & eventID==3)));


%%% bar plot of speed in different conditions (for grant)
figure
bar(spd)
hold on
errorbar([1 2 3],spd,spd_std,'ko','LineWidth',1)
set(gca,'xticklabels',{'CW','CCW','stop'})
ylabel('deg / sec')
xlim([0.25 3.75]); ylim([-45 45])
set(gca,'FontSize',12)


%%% create "wrapped" head angle, getting rid of discontinuities at +/-360deg
dth = diff(th);
dth(dth>90) = dth(dth>90)-360;
dth(dth<-90) = dth(dth<-90)+360;
th_wrap = cumsum(dth)
figure
plot((1:length(th_wrap))/30,th_wrap)
hold on
plot((1:length(stimV))/30,stimV);
xlim([0 120])

%%% speed and wrapped angle
stimV_short = stimV;
t= (1:length(stimV_short))/30-25;
figure
plot(t(stimV_short==50), th_wrap(stimV_short==50),'r.');
hold on
plot(t(stimV_short==-50), th_wrap(stimV_short==-50),'b.');
xlim([0 90])
legend('cw','ccw')

%%% adding lines for stim changes
t= (1:length(th_wrap))/30-25;
figure
plot(t,th_wrap)
hold on
for i = 1:length(cwstartF);
    plot([t(cwstartF(i)) t(cwstartF(i))],[-300 300],'r');

end
for i = 1:length(ccwstartF);
    plot([t(ccwstartF(i)) t(ccwstartF(i))],[-300 300],'b');    
end
xlim([0 90])
legend('th','cw','ccw')


%%% adding colored patches (final figure for grant)
t0 = 25;  %%% start time for plot (for grant, started at t=25sec to remove some short stim changes
t= (1:length(th_wrap))/30-t0; 
figure

hold on
for i = 1:length(eventF)-1;
    if eventID(i)==1
        col = 'b';
    elseif eventID(i)==2;
        col = 'g';
    elseif eventID(i)==3;
        col = 'r';
    end
    patch([t(eventF(i)) t(eventF(i+1)) t(eventF(i+1)) t(eventF(i))],[500 500 -500 -500],col,'FaceAlpha',0.25, 'EdgeColor','none')
end
plot(t,th_wrap,'k','lineWidth',2)
xlabel('secs'); ylabel('head angle (deg)')
xlim([0 88]); ylim([-250 400])
