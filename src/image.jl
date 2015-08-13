import Base: size, ndims, Array, length
using AffineTransforms
using Grid


type NeuroImage{T, P, N} <: AbstractArray{T, N}
    values :: AbstractArray{T, N}
    coordsys :: CoordSystem
    transformations :: CoordmapSet
    function NeuroImage(values :: AbstractArray, coordsys :: CoordSystem, transformations :: CoordmapSet) 
        if P > N
            error("Number of physical dimensions $P is too high for $N-dimensional array")
        end
        new(values, coordsys, transformations)
    end
end
NeuroImage{T, N}(values :: AbstractArray{T, N}, phys_ndims :: Int, coordsys :: CoordSystem, transformations :: CoordmapSet) = NeuroImage{T, phys_ndims, N}(values, coordsys, transformations)
NeuroImage(values :: AbstractArray, phys_ndims :: Int, coordsys :: CoordSystem) = NeuroImage(values, phys_ndims, coordsys, CoordmapSet())
NeuroImage(values, phys_ndims :: Int) = NeuroImage(values, phys_ndims, voxel())
NeuroImage(values :: AbstractArray) = NeuroImage(values, min(3, ndims(values)))
Array(img :: NeuroImage) = img.values
coordsys(img :: NeuroImage) = img.coordsys

size(img :: NeuroImage) = size(img.values)
size(img :: NeuroImage, d) = size(img.values, d)
ndims{T, P, N}(img :: NeuroImage{T, P, N}) = N
length(img :: NeuroImage) = length(img.values)
endof(img :: NeuroImage) = endof(img.values)

phys_ndims{T, P, N}(img :: NeuroImage{T, P, N}) = P
phys_size(img :: NeuroImage) = size(img)[1:phys_ndims(img)]

typealias RangeIndex Union(Int, Range{Int}, UnitRange{Int}, Colon)
function getindex{T, P, N}(img :: NeuroImage{T, P, N}, inds :: RangeIndex...)
    new_arr = img.values[inds...]
    if all(dim -> dim > ndims(new_arr) ? true : size(new_arr, dim) == 1, 1:phys_ndims(img))
        return new_arr
    end
    if ndims(new_arr) < P
        new_arr = reshape(new_arr, ntuple(P, ix -> ix <= ndims(new_arr) ? size(new_arr, ix) : 1))
    end
    new_coordsys = voxel()
    new_trans = union(img.transformations, CoordmapSet([Coordmap(new_coordsys, img.coordsys, inds_to_transform(phys_ndims(img), inds[1:phys_ndims(img)]...))]))
    NeuroImage{T, P, ndims(new_arr)}(new_arr, new_coordsys, new_trans)
end

setindex!(img :: NeuroImage, values, inds :: RangeIndex...) = img.values[inds...] = values

find_transformation(from_img :: NeuroImage, target_coordsys :: CoordSystem) = find_transformation(from_img.transformations, from_img.coordsys, target_coordsys, phys_ndims(from_img))
find_transformation(from_coordsys :: CoordSystem, target_img :: NeuroImage) = find_transformation(target_img.transformations, from_coordsys, target_img.coordsys, phys_ndims(target_img))
find_transformation(from_img :: NeuroImage, target_img :: NeuroImage) = find_transformation(union(from_img.transformations, target_img.transformations), from_img.coordsys, target_img._coordsys, phys_ndims(from_img))

location(img :: NeuroImage, coordsys :: CoordSystem, inds :: Integer...) = find_transformation(img, coordsys) * [dim > size(inds, 1) ? 1 : inds[dim] for dim=1:phys_ndims(img)]
location(img :: NeuroImage, inds :: Integer...) = [dim > size(inds, 1) ? 1 : inds[dim] for dim=1:phys_ndims(img)]

function physical_interpolator{T, P, N}(img :: NeuroImage{T, P, N}, interpolator=InterpLinear, boundary_condition=BCnil)
    non_phys_size = ntuple(N - P, dim -> size(img, dim + P))
    arr = reshape(Array(img), ntuple(P + 1, dim -> dim > P ? prod(non_phys_size) : size(img, dim)))
    grid_size = size(arr)
    grid_strides = strides(arr)
    strd = stride(arr, P+1)
    interp_invert!(arr, boundary_condition, interpolator, 1:P)   
    ic = InterpGridCoefs(eltype(arr), interpolator, grid_size[1:P], grid_strides[1:P])    
    function return_arr(index :: Vector)
        set_position(ic, boundary_condition, false, false, index)
        spec = Array(Float64, non_phys_size)
        index = 1
        for i = 1:prod(non_phys_size)
            spec[i] = interp(ic, arr, index)
            index += strd
        end
        spec
    end
    return_arr
end

function tuple_index{T <: Integer, S, N}(A :: AbstractArray{S, N}, index :: T)
    Astrides = strides(A)
    res = Array(T, N)
    for dim in sortperm(Astrides)[end:-1,1]
        res[dim] = 1 + ifloor((index - 1) / Astrides[dim])
        index = index - res[dim] * Astrides[dim]
    end
    res
end

function ndindex{N}(sz :: NTuple{N, Integer})
    if N == 0
        produce(())
    else
        for subset in @task ndindex(sz[2:end])
            for index in 1:sz[1]
                produce(tuple(index, subset...))
            end
        end
    end
end

function ndindex{T, P, N}(img :: NeuroImage{T, P, N})
    addition = ntuple(N - P, dim -> :)
    for index in @task ndindex{P}(phys_size(image))
        produce(tuple(index..., addition...))
    end
end
        
function transform_to{T, P, N}(from_img :: NeuroImage{T, P, N}, target_coordsys :: CoordSystem, phys_size :: NTuple{P, Integer}, transformations :: CoordmapSet=CoordmapSet([]), interpolator :: Grid.DataType=InterpLinear)
    sz = tuple(phys_size..., size(from_img)[P+1:N]...)
    new_img = NeuroImage(Array(T, sz), P, target_coordsys, transformations)
    trans = find_transformation(union(from_img.transformations, transformations), target_coordsys, from_img.coordsys, P)
    func = physical_interpolator(from_img, interpolator, BCnan)
    nvalues = prod(sz[P+1:N])
    for index in @task ndindex(sz[1:P])
        new_img[index..., 1:end] = func(trans * [index...])
    end
    new_img
end

transform_to{T, S, P, N, M}(from_img :: NeuroImage{T, P, N}, target_img :: NeuroImage{S, P, M}) = transform_to(from_img, target_img.coordsys, phys_size(target_img), target_img.transformations)

push_transformation!(img :: NeuroImage, transformation, to_coordsys :: CoordSystem, from_coordsys :: CoordSystem) = push!(img.transformations, Coordmap(from_coordsys, to_coordsys, transformation))
push_transformation!(img :: NeuroImage, transformation, to_coordsys :: CoordSystem) = push_transformation!(img, transformation, to_coordsys, img.coordsys)
