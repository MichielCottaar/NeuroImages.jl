using NeuroImages
using Base.Test
using AffineTransforms

arr = rand(8, 6, 5, 4)
img = NeuroImage(arr)

# reference space tests
@test voxel() != voxel()
@test mni() == mni()
@test coordsys(:MNI) == mni()
@test coordsys(:Talairach) == talairach()
@test mni() != talairach()
@test coordsys(img) == coordsys(img)
@test coordsys(NeuroImage(arr)) != coordsys(img)

# basic subsampling and physical, non-physical dimension tracking
@test ndims(img) == 4
@test phys_ndims(img) == 3
@test phys_size(img) == (8, 6, 5)
@test size(img) == (8, 6, 5, 4)
@test phys_ndims(img[2, 4:6, :, :]) == 3
@test ndims(img[2, 4:6, :, :]) == 4
@test phys_ndims(img[:, :, :, 4]) == 3
@test ndims(img[:, :, : 4]) == 3
@test phys_ndims(img[:, :, 2, 3]) == 3
@test ndims(img[:, :, 2, 3]) == 3

@test Array(img[2, 4:6, :, :]) == arr[2, 4:6, :, :]

# Voxel behaviour
@test img[2, 3, 5, 1] == arr[2, 3, 5, 1]
@test img[2, 3, 5, 1:3] == arr[2, 3, 5, 1:3]
@test location(img, 2, 3, 5, 1) == [2, 3, 5]
@test location(img, 2, 3, 5) == [2, 3, 5]
@test location(img, 2, 3) == [2, 3, 1]
@test location(img, img.coordsys, 2, 3, 5) == [2, 3, 5]
sub_img = img[4, 3:5, :, :]
@test location(sub_img, img.coordsys, 1, 1, 2) == [4, 3, 2]
@test sub_img[1, 1, 2] == arr[4, 3, 2]
@test location(sub_img, img.coordsys, 1, 2, 1) == [4, 4, 1]
@test sub_img[1, 2, 1] == arr[4, 4, 1]
@test location(sub_img, sub_img.coordsys, 1, 2, 1) == [1, 2, 1]

# Transformations
sub_img = img[4, 3:5, :, :]
@test transform_to(sub_img, sub_img) == sub_img
@test transform_to(img, sub_img) == sub_img
@test transform_to(img[3:end, :, :, :], sub_img) == sub_img
@test transform_to(img[3:end, :, :, 1], sub_img) == sub_img[:, :, :, 1]
sub_img2 = transform_to(img[:, :, 2:end, :], sub_img)
@test sub_img[:, :, 2:end, :] == sub_img2[:, :, 2:end, :]
@test all(isnan(sub_img2[:, :, 1:1, :, :])) # this slice is out of range

# Multiple transformations including standard space
push_transformation!(img, tformscale(2., 3), mni())
img2 = NeuroImage(rand(7, 7, 2))
push_transformation!(img2, tformtranslate([2, 0, 3]), mni(), img2.coordsys)
as_img2 = transform_to(img, img2)
@test as_img2[2, 2, 1] == img[2, 1, 2]
@test as_img2[4, 4, 1] == img[3, 2, 2]
img3 = NeuroImage(Array(img2))
push_transformation!(img, tformtranslate([-2, 0, -3]), img3.coordsys, mni())
@test all((Array(as_img2) .== Array(transform_to(img, img3))) | isnan(Array(as_img2)))
as_img = transform_to(img2, img)
@test all((Array(as_img) .== Array(transform_to(img3, img))) | isnan(Array(as_img)))

@test location(img, mni(), 2, 3) == [4, 6, 2]
@test location(img, img3.coordsys, 2, 3) == [2, 6, -1]

img_unz = read_nifti(joinpath(dirname(@__FILE__), "sample_image.nii"))
img_gz = read_nifti(joinpath(dirname(@__FILE__), "sample_image.nii.gz"))

@test Array(img_unz) == Array(img_gz)
@test_approx_eq(img_unz[13, 14, 16], 632.32965)
@test_approx_eq(location(img_unz, mni(), 13, 14, 16), [48.75, -47.25,  40.5])
