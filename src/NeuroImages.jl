module NeuroImages

export NeuroImage, phys_ndims, phys_size, location, transform_to, push_transformation!, mni, talairach, voxel, scanner, coordsys, read_nifti

include("reference.jl")
include("transformation.jl")
include("image.jl")
include("io.jl")

end # module
