# Makefile to create supergrid.nc and interpolated_topog.nc
# To use:
#   module load python
#   setenv PYTHONPATH $cwd/MIDAS
#
# then
#   make all
# or
#   make supergrid.nc
#   make interpolated_topog.nc

all: supergrid.nc interpolated_topog.nc

# Grids
supergrid.nc: mercator_supergrid.nc ncap_supergrid.nc antarctic_spherical_supergrid.nc scap_supergrid.nc
	python merge_grids.py

mercator_supergrid.nc ncap_supergrid.nc antarctic_spherical_supergrid.nc scap_supergrid.nc:
	python create_grids.py

# Topography
interpolated_topog.nc: ncap_topog.nc ncap_topog_gebco.nc mercator_topog_gebco.nc scap_topog_bedmap2.nc so_topog_bedmap2.nc scap_topog_gebco.nc so_topog_gebco.nc
	python merge_topog_tiles.py

ncap_topog.nc: mercator_supergrid.nc ncap_supergrid.nc antarctic_spherical_supergrid.nc scap_supergrid.nc
	python create_topo.py ncap

ncap_topog_gebco.nc: mercator_supergrid.nc ncap_supergrid.nc antarctic_spherical_supergrid.nc scap_supergrid.nc
	python create_topo.py ncap --use_gebco

mercator_topog_gebco.nc: mercator_supergrid.nc ncap_supergrid.nc antarctic_spherical_supergrid.nc scap_supergrid.nc
	python create_topo.py mercator

scap_topog_bedmap2.nc so_topog_bedmap2.nc: mercator_supergrid.nc ncap_supergrid.nc antarctic_spherical_supergrid.nc scap_supergrid.nc
	python create_topo.py scap

scap_topog_gebco.nc so_topog_gebco.nc: mercator_supergrid.nc ncap_supergrid.nc antarctic_spherical_supergrid.nc scap_supergrid.nc
	python create_topo.py scap --use_gebco

MIDAS:
	git clone https://github.com/mjharriso/MIDAS.git

MIDAS/README.md:
	cd ../../../; git submodule update

md5sums.txt: supergrid.nc interpolated_topog.nc
	echo Grids > $@
	md5sum *supergrid.nc >> $@
	echo >> $@
	echo Topography >> $@
	md5sum *topog*.nc >> $@
