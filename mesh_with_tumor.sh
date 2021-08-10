#!/bin/bash
# Glioma Brain Fluid Flow
# This script creates the .h5 file (mesh) from the result of the recon_all command
# The tumor is marked with label 999 in the .mgz parcellation files


# Give the path to the subjects directory 
SUBJECTS_DIR=~/Documents

# Name of the subject
NAME=TC2_Intrathecal

# Create the tumor surface from the wmparc file that contains the marking of the tumor
PARC_FILE=wmparc_with_tumor_marks.mgz
mri_binarize --i $SUBJECTS_DIR/$NAME/mri/$PARC_FILE --match 999 --surf-smooth 2 --surf tumor-surface.stl
# create the surface of the ventricles
mri_binarize --i $SUBJECTS_DIR/$NAME/mri/$PARC_FILE --ventricles --surf-smooth 2 --surf ventricles.stl

# copy the surfaces of gray and white matter to the correct directory
cp $SUBJECTS_DIR/$NAME/surf/lh.white .
cp $SUBJECTS_DIR/$NAME/surf/rh.white .
cp $SUBJECTS_DIR/$NAME/surf/lh.pial .
cp $SUBJECTS_DIR/$NAME/surf/rh.pial .

# convert surfaces to .stl files
mris_convert ./lh.white white.stl
mris_convert ./rh.white white.stl
mris_convert ./lh.pial pial.stl
mris_convert ./rh.pial pial.stl

# improve the surfaces for the pial membrane and the white matter (filling holes and separating narrow spaces)
# if needed?
# improve the surface of the ventricles (close the ventricles)

# Meshing of the brain and creation of the subdomains
MESH_FILE=($NAME".mesh")
RESOLUTION=32
python3 merging_hemisphere_6_surfaces.py --meshfile $MESH_FILE --resolution $RESOLUTION
H5_FILE=($NAME".h5")
# convert it to mesh for fenics
python3 convert_to_dolfin_mesh.py --meshfile $MESH_FILE --hdf5file $H5_FILE