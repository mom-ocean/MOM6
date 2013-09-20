ZLIM=[-250 0];
interpcolormap('bgr')

ncload('prog.nc'); ncload('visc.nc');
% Check that the state is horizontally uniform
dt = abs(temp(:,:,2,1) - temp(:,:,1,1));
dt = dt + abs(temp(:,:,1,2) - temp(:,:,1,1));
dt = dt + abs(temp(:,:,2,2) - temp(:,:,1,1));
if max(dt(:)) > 0
 stats(dt);error
end
% Throw away the i,j indics since the state is horizontally uniform
allvars=whos;
for j = 1:length(allvars)
 q=allvars(j);
 if length(q.size)==3 && prod(q.size(2:3))==4
   eval( sprintf('%s = %s(:,1);',q.name,q.name) )
 end
 if length(q.size)==4 && prod(q.size(3:4))==4
   eval( sprintf('%s = %s(:,:,1);',q.name,q.name) )
 end
end

subplot(821)
plot(Time,KPP_uStar); xlim([0 max(Time)])
%xlabel('Time (days)')
ylabel('u* (m/s)')

subplot(822)
plot(Time,KPP_buoyFlux); xlim([0 max(Time)])
%xlabel('Time (days)')
ylabel('B (m^2/s^3)')

subplot(823)
plot(Time,temp(:,1)); xlim([0 max(Time)])
ylim([18 27])
%xlabel('Time (days)')
ylabel('SST (^oC)')

subplot(824)
plot(Time,salt(:,1)); xlim([0 max(Time)])
ylim([36.25 36.75])
%xlabel('Time (days)')
ylabel('SSS (ppt)')

subplot(423)
gcolor(temp',e',Time);ylim(ZLIM);colorbar
caxis([18 27])
%xlabel('Time (days)')
ylabel('z (m)')
title('\theta (^oC)')
hold on;plot(Time,-h_ML,'w');hold off

subplot(424)
gcolor(salt',e',Time);ylim(ZLIM);colorbar
caxis([36.25 36.75])
%xlabel('Time (days)')
ylabel('z (m)')
title('S (ppt)')
hold on;plot(Time,-h_ML,'w');hold off

subplot(425)
gcolor(KPP_dTdt',e',Time);ylim(ZLIM);colorbar
caxis([-1 1]*1e-5)
%xlabel('Time (days)')
ylabel('z (m)')
title('KPP dTdt non-local (K/s)')
hold on;plot(Time,-h_ML,'w');hold off

subplot(426)
gcolor(KPP_dSdt',e',Time);ylim(ZLIM);colorbar
caxis([-1 1]*5e-7)
%xlabel('Time (days)')
ylabel('z (m)')
title('KPP dSdt non-local (ppt/s)')
hold on;plot(Time,-h_ML,'w');hold off

subplot(427)
gcolor(KPP_Ksalt',e',Time);ylim(ZLIM);colorbar
stats(KPP_Ksalt,'KPP_Ksalt')
caxis([0 2e-2])
xlabel('Time (days)')
ylabel('z (m)')
title('KPP \kappa_s (m^2/s)')
hold on;plot(Time,-KPP_OBLdepth ,'w');hold off

subplot(428)
gcolor(KPP_N',e',Time);ylim(ZLIM);colorbar
stats(KPP_N,'KPP N')
caxis([0 2e-2])
xlabel('Time (days)')
ylabel('z (m)')
title('KPP N (1/s)')
hold on;plot(Time,-h_ML,'w');hold off
