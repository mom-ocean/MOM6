OUT=ocean_hgrid.nc sponge.nc temp_salt_ic.nc topog.nc

$(OUT):
	python build_grid.py
	python build_data.py

clean:
	rm -rf $(OUT)
