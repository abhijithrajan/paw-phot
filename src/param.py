#! /usr/bin/env python


pixel_scale = 		[0.288]
datadir = 		["data_binned"]
variable_aperture = 	["yes"]
airmass_detrend = 	["yes"]

cut_offs = [2.0,0.0,15,5,10.,10000.]	# [Cut off for each bin, MAD factor ref light curve cut off, initial ref sample number, min remaining refs]