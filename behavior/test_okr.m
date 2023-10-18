clear all

d = readtable('PXL_20231003_235458873DLC_resnet50_squid_okr_100223Oct5shuffle1_100000.csv');
p_thresh = 0.95
n_nan = 30;

dt = 1/60;

p = d{:,4:3:end};

figure
plot(p)
figure
for i = 1:11
    subplot(4,3,i);
    plot(p(:,i))
end


head(:,1:2) = d{:,20:21};
head_p = d{:,22};

figure
plot(head_p)
head(head_p<p_thresh,:) = NaN;


eyeL(:,1:2) = d{:,5:6};
eyeL_p = d{:,7};
figure
plot(eyeL_p)
eyeL(eyeL_p<p_thresh,:) = NaN;
for i = 1:2
    eyeL(:,i) = interpNan(eyeL(:,i),n_nan,'pchip');
end

eyeR(:,1:2) = d{:,14:15};
eyeR_p = d{:,16};
figure
plot(eyeR_p)
eyeR(eyeR_p<p_thresh,:) = NaN;
for i = 1:2
    eyeR(:,i) = interpNan(eyeR(:,i),n_nan,'pchip');
end



eye_mid = 0.5*(eyeL+eyeR);
sum(~isnan(eye_mid(:,1)))/length(eye_mid(:,1))


th = atan2(eyeR(:,2)-eyeL(:,2),eyeR(:,1)-eyeL(:,1));
t = dt*(1:length(th));
figure
plot(t,th,'.')
xlim([0 60])
%xlim(dt*[0 10^4])

offset = 30;
dth = th((1+offset):end) - th(1:(end-offset));
dth(dth>pi) = dth(dth>pi)-2*pi;
dth(dth<-pi) = dth(dth<-pi)+2*pi;
figure
t = dt*(1:length(dth));
plot(t,dth); hold on
plot([1 t(end)],[0 0]);
xlim([0 60]); ylim([-1.5 1.5])

trange = 1:1000;

figure
plot(head(trange,1),head(trange,2))
hold on
%plot(eyeL(trange,1),eyeL(trange,2))
%plot(eyeR(trange,1),eyeR(trange,2))
plot(eye_mid(trange,1),eye_mid(trange,2))
axis equal