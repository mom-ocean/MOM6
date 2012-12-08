nc=netcdf('MOM_IC.nc');
x=nc{'lonh'}(:);
eta=squeeze( nc{'eta'}(:,1,:) );
salt=squeeze( nc{'Salt'}(1,:,1,:) );
gcolor(salt,eta,x); colorbar
xlabel(sprintf('X (%s)',nc{'lonh'}.units(:)))
ylabel(sprintf('Z (%s)',nc{'eta'}.units(:)))
title(sprintf('%s (%s)',nc{'Salt'}.long_name(:),nc{'Salt'}.units(:)))
close(nc)
print -dpng Fig_a.png

nc=netcdf('prog.nc');
eta=squeeze( nc{'e'}(10,:,1,:) );
salt=squeeze( nc{'salt'}(10,:,1,:) );
gcolor(salt,eta,x); colorbar
xlabel(sprintf('X (%s)',nc{'xh'}.units(:)))
ylabel(sprintf('Z (%s)',nc{'e'}.units(:)))
title(sprintf('%s (%s)',nc{'salt'}.long_name(:),nc{'salt'}.units(:)))
close(nc)
print -dpng Fig_b.png
