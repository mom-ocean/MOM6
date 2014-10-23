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

SHELL=tcsh -f

all: input forcing obs
	md5sum -c md5sums.txt

showenv:
	env
	-set
	-module list
	which python
	-python --version

input: ocean_hgrid.nc ocean_topog.nc basin_codes.nc

forcing: salt_restore.nc seawifs_1998-2006_smoothed_2X.nc tidal_amplitude.nc

obs: WOA05_ptemp_salt_monthly.nc WOA05_ptemp_salt_annual.nc

# Grids
supergrid.nc: mercator_supergrid.nc ncap_supergrid.nc antarctic_spherical_supergrid.nc scap_supergrid.nc local
	unlimit stacksize; setenv PYTHONPATH ./local/lib/python; python merge_grids.py

mercator_supergrid.nc ncap_supergrid.nc antarctic_spherical_supergrid.nc scap_supergrid.nc: local
	unlimit stacksize; setenv PYTHONPATH ./local/lib/python; python create_grids.py

# Topography
interpolated_topog.nc: ncap_topog.nc ncap_topog_gebco.nc mercator_topog_gebco.nc scap_topog_bedmap2.nc so_topog_bedmap2.nc scap_topog_gebco.nc so_topog_gebco.nc local
	unlimit stacksize; setenv PYTHONPATH ./local/lib/python; python merge_topog_tiles.py

ncap_topog.nc: mercator_supergrid.nc ncap_supergrid.nc antarctic_spherical_supergrid.nc scap_supergrid.nc IBCAO_V3_500m_RR.grd local
	unlimit stacksize; setenv PYTHONPATH ./local/lib/python; python create_topo.py ncap

ncap_topog_gebco.nc: mercator_supergrid.nc ncap_supergrid.nc antarctic_spherical_supergrid.nc scap_supergrid.nc GEBCO_08_v1.nc local
	unlimit stacksize; setenv PYTHONPATH ./local/lib/python; python create_topo.py ncap --use_gebco

mercator_topog_gebco.nc: mercator_supergrid.nc ncap_supergrid.nc antarctic_spherical_supergrid.nc scap_supergrid.nc GEBCO_08_v1.nc local
	unlimit stacksize; setenv PYTHONPATH ./local/lib/python; python create_topo.py mercator

scap_topog_bedmap2.nc so_topog_bedmap2.nc: mercator_supergrid.nc ncap_supergrid.nc antarctic_spherical_supergrid.nc scap_supergrid.nc bedmap2.nc local
	unlimit stacksize; setenv PYTHONPATH ./local/lib/python; python create_topo.py scap

scap_topog_gebco.nc so_topog_gebco.nc: mercator_supergrid.nc ncap_supergrid.nc antarctic_spherical_supergrid.nc scap_supergrid.nc GEBCO_08_v1.nc local
	unlimit stacksize; setenv PYTHONPATH ./local/lib/python; python create_topo.py scap --use_gebco

# Sets char tile='tile1'
ocean_hgrid.nc: supergrid.nc
	ncks -h -d ny,80, -d nyp,80, supergrid.nc ocean_hgrid.nc
	./changeChar.py ocean_hgrid.nc tile tile1

topog.nc: edit_topog.nc interpolated_topog.nc
	./ice9.py edit_topog.nc --output topog.nc

ocean_topog.nc: topog.nc
	ncks -h -d ny,40,  topog.nc ocean_topog.nc
	ncap -h -s 'jEdit=jEdit-40' --overwrite ocean_topog.nc ocean_topog.nc
	./addDimension.py ocean_topog.nc ntiles 1
# Removes rows 1-40
# Adjusts jEdit values
# Adds back the dimension "ntiles"

basin_codes.nc: ocean_topog.nc
	unlimit stacksize; setenv PYTHONPATH ./local/lib/python; python make_basin_mask.py
	ncatted -h -a flag_meanings,basin,c,c,'1:Southern Ocean, 2:Atlantic Ocean, 3:Pacific Ocean, 4:Arctic Ocean, 5:Indian Ocean, 6:Mediterranean Sea, 7:Black Sea, 8:Hudson Bay, 9:Baltic Sea, 10:Red Sea, 11:Persian Gulf' basin_codes.nc
	ncatted -h -a flag_values,basin,c,c,'1,2,3,4,5,6,7,8,9,10,11' basin_codes.nc

WOA05_ptemp_salt_annual.nc: WOA05_ptemp_salt_monthly.nc
	ncra -h -O $< $@

WOA05_ptemp_salt_monthly.nc: ocean_hgrid.nc ocean_topog.nc /archive/gold/datasets/obs/WOA05_pottemp_salt.nc interpWOA05.py local
	unlimit stacksize; setenv PYTHONPATH ./local/lib/python; python interpWOA05.py
	ncatted -h -a modulo,TIME,c,c,' ' WOA05_ptemp_salt_monthly.nc

salt_restore.nc: ocean_hgrid.nc ocean_topog.nc /archive/gold/datasets/obs/WOA05_pottemp_salt.nc interpSaltRestore.py local
	unlimit stacksize; setenv PYTHONPATH ./local/lib/python; python interpSaltRestore.py
	ncatted -h -a modulo,TIME,c,c,' ' salt_restore.nc

seawifs_1998-2006_smoothed_2X.nc: ocean_hgrid.nc ocean_topog.nc /archive/gold/datasets/global/siena_201204/INPUT/seawifs_1998-2006_GOLD_smoothed_2X.nc interpCHL.py local
	unlimit stacksize; setenv PYTHONPATH ./local/lib/python; python interpCHL.py
	ncatted -h -a modulo,TIME,c,c,' ' seawifs_1998-2006_smoothed_2X.nc

tidal_amplitude.nc: DATA interpTides.py local
	unlimit stacksize; setenv PYTHONPATH ./local/lib/python; python interpTides.py

DATA: tpxo7_atlas_netcdf.tar.Z
	tar xzvf tpxo7_atlas_netcdf.tar.Z
	touch DATA

tpxo7_atlas_netcdf.tar.Z:
	wget ftp://ftp.oce.orst.edu/dist/tides/Global/tpxo7_atlas_netcdf.tar.Z

MIDAS:
	git clone https://github.com/mjharriso/MIDAS.git

local: MIDAS 
	-rm -rf $</build/*
	mkdir -p $@
	cd $<; make -f Makefile_GFDL INSTALL_PATH=../local
	touch $@

bedmap2.nc GEBCO_08_v1.nc IBCAO_V3_500m_RR.grd:
	cp /archive/gold/datasets/topography/$@ .

md5sums.txt: ocean_hgrid.nc antarctic_spherical_supergrid.nc mercator_supergrid.nc ncap_supergrid.nc scap_supergrid.nc supergrid.nc interpolated_topog.nc mercator_topog_gebco.nc ncap_topog_gebco.nc ncap_topog.nc scap_topog_bedmap2.nc scap_topog_gebco.nc so_topog_bedmap2.nc so_topog_gebco.nc ocean_topog.nc basin_codes.nc salt_restore.nc seawifs_1998-2006_smoothed_2X.nc tidal_amplitude.nc WOA05_ptemp_salt_monthly.nc WOA05_ptemp_salt_annual.nc
	echo Grids > $@
	md5sum *supergrid.nc ocean_hgrid.nc >> $@
	echo >> $@
	echo Topography >> $@
	md5sum *topog*.nc >> $@
	md5sum basin_codes.nc >> $@
	echo >> $@
	echo Data >> $@
	md5sum salt_restore.nc seawifs_1998-2006_smoothed_2X.nc tidal_amplitude.nc >> $@
	echo >> $@
	echo Obs >> $@
	md5sum WOA05_ptemp_salt_{monthly,annual}.nc >> $@