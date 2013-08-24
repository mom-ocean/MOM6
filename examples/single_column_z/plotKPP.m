ZLIM=[-250 0];
interpcolormap('bgr')

ncload('prog.nc'); ncload('visc.nc')
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

subplot(321)
gcolor(temp',e',Time);ylim(ZLIM);colorbar
caxis([18 21])
xlabel('Time (days)');ylabel('z (m)')
title('\theta (^oC)')
hold on;plot(Time,-KPP_OBLdepth,'w');hold off

subplot(322)
gcolor(salt',e',Time);ylim(ZLIM);colorbar
caxis([36 37])
xlabel('Time (days)');ylabel('z (m)')
title('S (ppt)')
hold on;plot(Time,-KPP_OBLdepth,'w');hold off

subplot(323)
gcolor(KPP_N',e',Time);ylim(ZLIM);colorbar
stats(KPP_N,'KPP_N')
caxis([0 2e-2])
xlabel('Time (days)');ylabel('z (m)')
title('KPP N (1/s)')
hold on;plot(Time,-KPP_OBLdepth,'w');hold off

subplot(324)
gcolor(Kd_interface',e',Time);ylim(ZLIM);colorbar
stats(Kd_interface,'Kd_interface')
caxis([0 2e-2])
xlabel('Time (days)');ylabel('z (m)')
title('Kd interface (m^2/s)')
hold on;plot(Time,-KPP_OBLdepth,'w');hold off

subplot(325)
gcolor(KPP_Ksalt',e',Time);ylim(ZLIM);colorbar
stats(KPP_Ksalt,'KPP_Ksalt')
caxis([0 2e-2])
xlabel('Time (days)');ylabel('z (m)')
title('KPP \kappa_s (m^2/s)')
hold on;plot(Time,-KPP_OBLdepth ,'w');hold off

subplot(6,2,10)
plot(Time,KPP_uStar)
xlabel('Time (days)');ylabel('u* (m/s)')

subplot(6,2,12)
plot(Time,KPP_buoyFlux)
xlabel('Time (days)');ylabel('B (m^2/s^3)')
