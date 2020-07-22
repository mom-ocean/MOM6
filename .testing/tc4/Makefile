ocean_hgrid.nc topog.nc temp_salt_ic.nc sponge.nc:
	python2 -m virtualenv local-env
	. local-env/bin/activate \
	  && pip install netCDF4 \
	  && python build_grid.py \
	  && python build_data.py
