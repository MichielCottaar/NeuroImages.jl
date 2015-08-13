using NIfTI: niread, isgz, getaffine

sform_code_values = [:ScannerAnatomical, :ScannerAligned, :Talairach, :MNI]
function sform_code(value :: Integer)
    if value == 0 || value > length(sform_code_values)
        coordsys(:Unknown)
    else
        coordsys(sform_code_values[value])
    end   
end

function read_nifti(filename)
    io = open(filename)
    niv = isgz(io) ? niread(filename) : niread(filename, mmap=true)
    nifti_space = coordsys(:ZeroBased)
    julia_space = coordsys(:OneBased)
    img = NeuroImage(niv.raw, 3, julia_space)
    affine = getaffine(niv.header)
    push_transformation!(img, tformtranslate([-1, -1, -1]), nifti_space, julia_space)
    push_transformation!(img, AffineTransform(affine[1:3, 1:3], affine[1:3, 4]), sform_code(niv.header.sform_code), nifti_space)
    img
end

    
