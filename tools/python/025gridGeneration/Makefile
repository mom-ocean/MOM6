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

all: md5sums.txt ocean_hgrid.nc

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
ocean_topog.nc: topog.nc
	ncks -h -d ny,40,  topog.nc ocean_topog.nc                             # Removes rows 1-40
	ncap -h -s 'jEdit=jEdit-40' --overwrite ocean_topog.nc ocean_topog.nc  # Adjusts jEdit values
	./addDimension.py ocean_topog.nc ntiles 1                              # Adds back the dimension "ntiles"

MIDAS:
	git clone https://github.com/mjharriso/MIDAS.git

MIDAS/README.md:
	cd ../../../; git submodule update

bedmap2.nc GEBCO_08_v1.nc IBCAO_V3_500m_RR.grd:
	ln -s /archive/gold/datasets/topography/$@ .

md5sums.txt: ocean_hgrid.nc antarctic_spherical_supergrid.nc mercator_supergrid.nc ncap_supergrid.nc scap_supergrid.nc supergrid.nc interpolated_topog.nc mercator_topog_gebco.nc ncap_topog_gebco.nc ncap_topog.nc scap_topog_bedmap2.nc scap_topog_gebco.nc so_topog_bedmap2.nc so_topog_gebco.nc
	echo Grids > $@
	md5sum *supergrid.nc ocean_hgrid.nc >> $@
	echo >> $@
	echo Topography >> $@
	md5sum *topog*.nc >> $@

