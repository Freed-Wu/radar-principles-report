#!/usr/bin/env zsh
dirs=(fig/beat fig/mti fig/doppler fig/mtd)
mkdir dirs
octave lst/main.m
for dir in $dirs; do
	for file in `ls $dir/*.tex`; do
		sed -i '2s/}/, ctex}/' $file
		sed -i 's='$dir'/==g' $file
		latexmk -cd -pvc- $file
		rm $file -f
	done
done
