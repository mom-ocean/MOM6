all: ocean_hgrid.nc topog.nc temp_salt_ic.nc sponge.nc

../local-env:
	python3 -m venv ../local-env
	source ../local-env/bin/activate \
	  && pip3 install wheel \
	  && pip3 install cython \
	  && pip3 install numpy \
	  && pip3 install netCDF4

ocean_hgrid.nc topog.nc: ../local-env
	source ../local-env/bin/activate \
	  && python build_grid.py

temp_salt_ic.nc sponge.nc: ../local-env
	source ../local-env/bin/activate \
	  && python build_data.py

clean:
	rm -rf local-env ocean_hgrid.nc topog.nc temp_salt_ic.nc sponge.nc
