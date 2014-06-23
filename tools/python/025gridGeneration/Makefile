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

SHELL=tcsh

all: md5sums.txt ocean_hgrid.nc ocean_topog.nc basin_mask.nc

forcing: salt_restore.nc seawifs_1998-2006_smoothed_2X.nc tidal_amplitude.nc

# Grids
supergrid.nc: mercator_supergrid.nc ncap_supergrid.nc antarctic_spherical_supergrid.nc scap_supergrid.nc
	python merge_grids.py

mercator_supergrid.nc ncap_supergrid.nc antarctic_spherical_supergrid.nc scap_supergrid.nc:
	python create_grids.py

# Topography
interpolated_topog.nc: ncap_topog.nc ncap_topog_gebco.nc mercator_topog_gebco.nc scap_topog_bedmap2.nc so_topog_bedmap2.nc scap_topog_gebco.nc so_topog_gebco.nc
	python merge_topog_tiles.py

ncap_topog.nc: mercator_supergrid.nc ncap_supergrid.nc antarctic_spherical_supergrid.nc scap_supergrid.nc IBCAO_V3_500m_RR.grd
	python create_topo.py ncap

ncap_topog_gebco.nc: mercator_supergrid.nc ncap_supergrid.nc antarctic_spherical_supergrid.nc scap_supergrid.nc GEBCO_08_v1.nc
	python create_topo.py ncap --use_gebco

mercator_topog_gebco.nc: mercator_supergrid.nc ncap_supergrid.nc antarctic_spherical_supergrid.nc scap_supergrid.nc GEBCO_08_v1.nc
	python create_topo.py mercator

scap_topog_bedmap2.nc so_topog_bedmap2.nc: mercator_supergrid.nc ncap_supergrid.nc antarctic_spherical_supergrid.nc scap_supergrid.nc bedmap2.nc
	python create_topo.py scap

scap_topog_gebco.nc so_topog_gebco.nc: mercator_supergrid.nc ncap_supergrid.nc antarctic_spherical_supergrid.nc scap_supergrid.nc GEBCO_08_v1.nc
	python create_topo.py scap --use_gebco

ocean_hgrid.nc: supergrid.nc
	ncks -h -d ny,80, -d nyp,80, supergrid.nc ocean_hgrid.nc
	./changeChar.py ocean_hgrid.nc tile tile1                              # Sets char tile='tile1'

topog.nc: edit_topog.nc
	./ice9.py edit_topog.nc --output topog.nc

ocean_topog.nc: topog.nc
	ncks -h -d ny,40,  topog.nc ocean_topog.nc
	ncap -h -s 'jEdit=jEdit-40' --overwrite ocean_topog.nc ocean_topog.nc
	./addDimension.py ocean_topog.nc ntiles 1
# Removes rows 1-40
# Adjusts jEdit values
# Adds back the dimension "ntiles"

basin_mask.nc: ocean_topog.nc
	./make_basin_mask.py

salt_restore.nc: local ocean_hgrid.nc ocean_topog.nc /archive/gold/datasets/obs/WOA05_pottemp_salt.nc
	unlimit stacksize; setenv PYTHONPATH ./local/lib/python; python interpSaltRestore.py
	module load nco; ncatted -h -a modulo,TIME,c,c,' ' salt_restore.nc

seawifs_1998-2006_smoothed_2X.nc: local ocean_hgrid.nc ocean_topog.nc /archive/gold/datasets/global/siena_201204/INPUT/seawifs_1998-2006_GOLD_smoothed_2X.nc
	unlimit stacksize; setenv PYTHONPATH ./local/lib/python; python interpCHL.py
	module load nco; ncatted -h -a modulo,TIME,c,c,' ' seawifs_1998-2006_smoothed_2X.nc

tidal_amplitude.nc: DATA
	unlimit stacksize; setenv PYTHONPATH ./local/lib/python; python interpTides.py

DATA: tpxo7_atlas_netcdf.tar.Z
	tar xzvf tpxo7_atlas_netcdf.tar.Z
	touch DATA

tpxo7_atlas_netcdf.tar.Z:
	wget ftp://ftp.oce.orst.edu/dist/tides/Global/tpxo7_atlas_netcdf.tar.Z

MIDAS:
	git clone https://github.com/mjharriso/MIDAS.git

local: MIDAS
	cd MIDAS/fms;tar xvf fms_tikal_201312.tar
	module load python netcdf/4.2 intel_compilers; cd MIDAS ; (cd fms_build;./build_fms.csh); python setup.py config_fc --f90flags="-i4 -r8 -DPY_SOLO" --fcompiler=intelem build
	cd MIDAS; python setup.py install --home=../local

MIDAS/README.md:
	cd ../../../; git submodule update

bedmap2.nc GEBCO_08_v1.nc IBCAO_V3_500m_RR.grd:
	ln -s /archive/gold/datasets/topography/$@ .

md5sums.txt: ocean_hgrid.nc antarctic_spherical_supergrid.nc mercator_supergrid.nc ncap_supergrid.nc scap_supergrid.nc supergrid.nc interpolated_topog.nc mercator_topog_gebco.nc ncap_topog_gebco.nc ncap_topog.nc scap_topog_bedmap2.nc scap_topog_gebco.nc so_topog_bedmap2.nc so_topog_gebco.nc ocean_topog.nc
	echo Grids > $@
	md5sum *supergrid.nc ocean_hgrid.nc >> $@
	echo >> $@
	echo Topography >> $@
	md5sum *topog*.nc >> $@

