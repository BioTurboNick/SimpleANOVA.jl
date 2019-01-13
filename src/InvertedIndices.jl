"""
    InvertedIndex(idx)
    Not(idx)

Construct an inverted index, selecting all indices not in the passed `idx`.

Upon indexing into an array, the `InvertedIndex` behaves like a 1-dimensional
collection of the indices of the array that are not in `idx`. Bounds are
checked to ensure that all indices in `idx` are within the bounds of the array
— even though they are skipped. The `InvertedIndex` behaves like a
1-dimensional collection of its inverted indices. If `idx` spans multiple
dimensions (like a multidimensional logical mask or `CartesianIndex`), then the
inverted index will similarly span multiple dimensions.

Copied and trimmed down from InvertedIndices.jl repository
"""
struct InvertedIndex{T}
    skip::T
end
const Not = InvertedIndex

import Base.trues
# This is a little tricky because trues isn't indices-aware, but we also don't
# want to use fill(true) in the 1-indexed case since we want to favor BitArrays
@inline trues(tup::Tuple{Vararg{Base.OneTo}}) = Base.trues(map(Base.unsafe_length, tup))
@inline trues(tup::Tuple{Vararg{Base.AbstractUnitRange}}) = fill(true, tup)

@inline function Base.to_indices(A, inds, I::Tuple{InvertedIndex, Vararg{Any}})
    v = trues(spanned_indices(inds, I))
    v[I[1].skip] = false
    to_indices(A, inds, (v, Base.tail(I)...))
end

# Determining the indices that the InvertedIndex spans is tricky due to partial
# linear indexing. Lean on `Base.uncolon` until the deprecation goes through.
@inline spanned_indices(inds, I::Tuple{InvertedIndex,Vararg{Any}}) = (Base.uncolon(inds, (:, Base.tail(I)...)).indices,)
