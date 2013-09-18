#! /usr/bin/env python
pixel_scale = 		[0.288]
rootdir = [data/bd/T-Y] # i.e "/home/user/data" Basically the path where the paw.py file is.
datadir = 		["data"]
variable_aperture = 	["yes"]
airmass_detrend = 	["yes"]

cut_offs = [2.0,0.0,15,5,10.,10000.]	# [Cut off for each bin, MAD factor ref light curve cut off, initial ref sample number, min remaining refs]

ignore_objects = [34,38,51,56,60,63,72,77]
