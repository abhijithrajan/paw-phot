#! /usr/bin/env python

import sys, shutil

import pyraf.iraf as iraf
from os  import chdir, path, makedirs


import os.path
from pyraf import iraf


def cd(target):
    chdir(target)
    iraf.cd(target)

def copy_to_dir(orig, target):
    filename = path.basename(orig)
    shutil.copyfile(orig, path.join(target, filename))

def move_to_dir(orig, target):
    filename = path.basename(orig)
    shutil.move(orig, path.join(target, filename))

def mkdir(dirname):
   "Creates a directory if it didn't exist before"
   if not os.path.isdir(dirname):
      os.mkdir(dirname)

def rmdir(dirname):
   "Remove a directory if it exist"
   if os.path.isdir(dirname):
      shutil.rmtree(dirname)

def bold(text):
   bold_code = "\033[1m"
   reset_code = "\033[0;0m"
   return bold_code + text + reset_code

def underline(text):
   underline_code = "\033[4m"
   reset_code = "\033[0;0m"
   return underline_code + text + reset_code

def blue(text):
   red = "\033[0;34m"
   reset_code = "\033[0;0m"
   return red + text + reset_code

def green(text):
   red = "\033[0;32m"
   reset_code = "\033[0;0m"
   return red + text + reset_code
   
def red(text):
   red = "\033[0;31m"
   reset_code = "\033[0;0m"
   return red + text + reset_code

def yellow(text):
   red = "\033[0;33m"
   reset_code = "\033[0;0m"
   return red + text + reset_code

def error(text):
   red = "\033[1;31m"
   reset_code = "\033[0;0m"
   return red + text + reset_code

def check_dir(directory):
   if not os.path.exists(directory):
     print error("\nError:\n")+"Cannot find the directory named: "+bold(directory)+"\n\nPlease run the sort function again.\n"
     sys.exit()

def copy_directory_if_exist(object_dir, fits_file, directory, colour):
   if os.path.exists(directory+"/"+colour):
     copy_to_dir(path.join(object_dir, fits_file), directory+"/"+colour)
     if colour == "blue":
       print blue(fits_file)+" ---> "+directory+"/"+blue(colour)+"\n"
     elif colour == "green":
       print green(fits_file)+" ---> "+directory+"/"+green(colour)+"\n"
     elif colour == "red":
       print red(fits_file)+" ---> "+directory+"/"+red(colour)+"\n"
     elif colour == "R500R":
       print yellow(fits_file)+" ---> "+directory+"/"+yellow(colour)+"\n"                     
          
def create_dir(directory):
   if not os.path.exists(directory):
     print "Creating \'"+directory+"\' directory\n"
     os.makedirs(directory)
     
def file_found(filetype,filename,colour):
  if filetype == "flat" and colour=="R500R":
    print "Flat file found: "+red(filename)+"\n"
  if filetype == "res_flat" and colour=="grisms":
    print "Flat file found: "+red(filename)+"\n"
  if filetype == "illum_flat" and colour=="R500B":
    print "Illum. flat file found: "+yellow(filename)+"\n"
  if filetype == "object" and colour=="R500B":
    print "Science files found.\n"
  
  
  
  
  
