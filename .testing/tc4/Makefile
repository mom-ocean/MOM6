# TODO: Resolve python path
# TODO: Identify if venv is available
PYTHON=python3
PIP=pip3

HAS_NUMPY=$(shell $(PYTHON) -c "import numpy" && echo "yes")
HAS_NETCDF4=$(shell $(PYTHON) -c "import netCDF4" && echo "yes")

OUT=ocean_hgrid.nc sponge.nc temp_salt_ic.nc topog.nc

all: $(OUT)

$(OUT):
ifneq ($(HAS_NUMPY), yes)
ifneq ($(HAS_NETCDF4), yes)
	$(PYTHON) -m venv local-env
	. local-env/bin/activate \
	  && $(PIP) install wheel \
	  && $(PIP) install cython \
	  && $(PIP) install numpy \
	  && $(PIP) install netCDF4
endif
endif
	$(PYTHON) build_grid.py
	$(PYTHON) build_data.py

clean:
	rm -rf local-env $(OUT)
