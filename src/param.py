#! /usr/bin/env python
pixel_scale = 		[0.288]
rootdir = 			["/data1/paw/bd/T-Y/binned"] # i.e "/home/user/data" Basically the path where the paw.py file is.
datadir = 			["data"]
image_name = 		["binned",".fits"]		# First and last parts of the .fits files.
variable_aperture = ["yes"]
airmass_detrend = 	["yes"]

cut_offs = [2.0,0.0,15,5,10.,10000.]	# [Cut off for each bin, MAD factor ref light curve cut off, initial ref sample number, min remaining refs]

ignore_objects = [34,38,51,56,60,63,72,77]
