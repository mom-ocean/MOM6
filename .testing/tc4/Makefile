# TODO: Identify if venv is available
HAS_NUMPY=$(shell python -c "import numpy" 2> /dev/null && echo "yes")
HAS_NETCDF4=$(shell python -c "import netCDF4" 2> /dev/null && echo "yes")

ifneq ($(HAS_NUMPY), yes)
ifneq ($(HAS_NETCDF4), yes)
USE_VENV=yes
endif
endif

OUT=ocean_hgrid.nc sponge.nc temp_salt_ic.nc topog.nc

all: $(OUT)

$(OUT):
ifeq ($(USE_VENV),yes)
	python3 -m venv local-env
	. local-env/bin/activate \
	  && pip3 install wheel \
	  && pip3 install cython \
	  && pip3 install numpy \
	  && pip3 install netCDF4 \
	  && python3 build_grid.py \
	  && python3 build_data.py
else
	python build_grid.py
	python build_data.py
endif

clean:
	rm -rf local-env $(OUT)
