# NeuroImages

[![Build Status](https://travis-ci.org/MichielCottaar/NeuroImages.jl.svg?branch=master)](https://travis-ci.org/MichielCottaar/NeuroImages.jl)

This repository defines an array-like NeuroImage type, which has been designed to represent a volumetric image of e.g. a brain. Currently the main addition to a regular array is that this image is assumed to lie in physical space with a well-defined coordinate system. The first N (typically 2 or 3) dimensions of the image are assumed to represent physical space with the remaining dimensions available for any other variations (e.g. time or RGB). The number of physical dimensions can be retreived by `phys_ndims(image)` or the sizes with `phys_size(image)`.

The current coordinate system of the image is stored in the coordsys. At present a coordinate system is represented by a name and whether it is absolute (e.g. MNI or Talairach produced from `mni()` or `talairach()`) or relative (e.g. voxel or scanner space produced from `voxel()` or `scanner()`). Any other coordinate systems can be produced by calling `coordsys`. If a subimage is created a new voxel coordinate system is created with an appropriate transformation to the total image.

Normally the images will be stored in voxel space, but an image can contain any number of transformations (stored in the `transformations` field) to other coordinate systems (currently only affine transformations are supported). When `transform_to` is called these transformations are automatically combined to obtain a master transformations, which can be used to convert between images in different coordinates system if possible. 

New images can be created simply by calling `NeuroImage(array)` or from a NIFTI file `read_nifti("filename")`.
