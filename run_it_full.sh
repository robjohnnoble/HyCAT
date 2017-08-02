#!/bin/bash

shopt -s extglob

for conc in 0 5 10 20 50 # seeds
	do
for i in {100..999} # seeds
	do

	outputdir="output_paper4_conc${conc}_seed${i}_withgrids"
	subdir_type="animation_frames_type"
	subdir_birth="animation_frames_birth"
	subdir_birth_type1="animation_frames_birth_type0"
	subdir_birth_type2="animation_frames_birth_type1"
	subdir_death="animation_frames_death"
	subdir_death_type1="animation_frames_death_type0"
	subdir_death_type2="animation_frames_death_type1"
	subdir_maxbirth="animation_frames_maxbirth"
	subdir_stem="animation_frames_stem"
	matsubdir="matrices"
	prosubdir="max_prolif_rates"

	cd "Results"
	mkdir $outputdir
	cd ..
	cp "parameters_paper${conc}.cfg" "Results/${outputdir}/parameters.cfg"
	cp HyCAT.c "Results/${outputdir}/HyCAT.c"
	cd "Results"
	cd $outputdir

	sed -i.bak "s/idum=0/idum=$i/g" parameters.cfg

	mkdir $subdir_type
	mkdir $subdir_birth
	mkdir $subdir_birth_type1
	mkdir $subdir_birth_type2
	mkdir $subdir_death
	mkdir $subdir_death_type1
	mkdir $subdir_death_type2
	mkdir $subdir_maxbirth
	mkdir $subdir_stem
	mkdir $matsubdir
	mkdir $prosubdir

	clang -O2 $(pkg-config --cflags libconfig) HyCAT.c -o HyCAT $(pkg-config --libs libconfig)

	chmod +x HyCAT
	./HyCAT

	convert -delay 1 -loop 0 "${subdir_type}/ani*.png" animation.gif

	convert -delay 1 -loop 0 "${subdir_birth}/occupied*.png" occupied.gif

	convert -delay 1 -loop 0 "${subdir_stem}/stem*.png" stem.gif

	convert -delay 1 -loop 0 "${subdir_maxbirth}/maxprolif*.png" maxprolif.gif

	gifs="$(find $subdir_type -maxdepth 1 -name 'ani*.png'|sort|xargs)"

	_CONVERT="convert "
	for gif in $gifs; do delay=$(echo ${gif##*/} | cut -d\_ -f3 | sed "s,.png,,"); _CONVERT="${_CONVERT} -delay $delay $gif "; done;
	# option -f specifies which field to extract
	# option -d specifies the field delimiter
	# sed "s,.png,," removes the extension ".png"

	_CONVERT="${_CONVERT} -loop 0 animation_variable.gif"
	#echo "Convert cmd: $_CONVERT"

	eval $_CONVERT

	cd ..
	cd ..
done
done