#!/bin/bash

shopt -s extglob

subdir_type="animation_frames_type"

gifs="$(find $subdir_type -maxdepth 1 -name 'ani*.png'|sort|xargs)"

_CONVERT="convert "
for gif in $gifs; do delay=$(echo ${gif##*/} | cut -d\_ -f3 | sed "s,.png,,"); _CONVERT="${_CONVERT} -delay $(($delay*10)) $gif "; done;
# option -f specifies which field to extract
# option -d specifies the field delimiter
# sed "s,.png,," removes the extension ".png"

_CONVERT="${_CONVERT} -loop 0 animation_variable_longer_gap.gif"
#echo "Convert cmd: $_CONVERT"

eval $_CONVERT