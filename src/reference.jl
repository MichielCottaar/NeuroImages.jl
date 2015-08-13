import Base.ndims, Base.inv
using AffineTransforms

abstract CoordSystem

immutable AbsoluteCoordSystem <: CoordSystem
    name :: Symbol
end

type RelativeCoordSystem <: CoordSystem
    name :: Symbol
end

absolute_coord = Set([:MNI, :Talairach])

coordsys(name :: Symbol, absolute :: Bool) = absolute ? AbsoluteCoordSystem(name) : RelativeCoordSystem(name)
coordsys(name :: Symbol) = coordsys(name, name in absolute_coord)
    
mni() = coordsys(:MNI, true)
talairach() = coordsys(:Talairach, true)
voxel() = coordsys(:voxel, false)
scanner() = coordsys(:scanner, false)


immutable Coordmap
    domain :: CoordSystem
    range :: CoordSystem
    transformation
end

function *(coord1 :: Coordmap, coord2 :: Coordmap)
    if coord1.domain != coord2.range
        error("Coordinate maps can not be composed because of inconsistent input, output coordinate systems")
    end
    Coordmap(coord2.domain, coord1.range, coord1.transformation * coord2.transformation)
end

inv(coord :: Coordmap) = Coordmap(coord.range, coord.domain, inv(coord.transformation))

ndims(coord :: Coordmap) = ndims(transformation)

typealias CoordmapSet Set{Coordmap}

function find_transformation(coord_set :: CoordmapSet, initial_coordsys :: CoordSystem, target_coordsys :: CoordSystem, ndims :: Integer)
    found_coordsys = Dict{CoordSystem, Vector{Coordmap}}([initial_coordsys], {[Coordmap(initial_coordsys, initial_coordsys, tformeye(ndims))]})
    previous_run = 0
    while found_coordsys.count > previous_run
        previous_run = found_coordsys.count
        for coordmap in coord_set
            if (coordmap.domain in keys(found_coordsys)) && ~(coordmap.range in keys(found_coordsys))
                found_coordsys[coordmap.range] = [[coordmap], found_coordsys[coordmap.domain]]
            elseif coordmap.range in keys(found_coordsys) && ~(coordmap.domain in keys(found_coordsys))
                found_coordsys[coordmap.domain] = [[inv(coordmap)], found_coordsys[coordmap.range]]
            end
            if target_coordsys in keys(found_coordsys)
                return prod(found_coordsys[target_coordsys]).transformation
            end
        end
        if target_coordsys in keys(found_coordsys)
            return prod(found_coordsys[target_coordsys]).transformation
        end
    end 
    error("No transformation to target coordinate system has been found")
end
