#!/bin/bash
# Glioma Brain Fluid Flow
# This script creates the .h5 file (mesh) from the result of the recon_all command
# The tumor is marked with label 999 in the .mgz parcellation files


# Give the path to the subjects directory 
SUBJECTS_DIR = ""
# Name of the subject
NAME = ""

# Create the tumor surface from the wmparc file that contains the marking of the tumor
PARC_FILE = ""
mri_binarize --i $SUBJECTS_DIR/mir/$PARC_FILE --match 999 --surf-smooth 2 --surf tumor-surface.stl
# create the surface of the ventricles
mri_binarize --i $SUBJECTS_DIR/mir/$PARC_FILE --ventricles --surf-smooth 2 --surf ventricles.stl

# improve the surfaces for the pial membrane and the white matter (filling holes and separating narrow spaces)

# improve the surface of the ventricles (close the ventricles)

# Meshing of the brain and creation of the subdomains

