abstract type Abstract_State end


#=====================================================================================================================================#
#======================================================== Single-Particle State ======================================================#
#=====================================================================================================================================#
"""
General Single-Particle State `Single_Particle_State{D}` <: Abstract_State`
---
with `D<:Int` the number of physical dof.
- Fields:
    - `dof_indices::Tuple`: tuple of indices to each physical dof. The order of each dof should be managed by the user
"""
struct Single_Particle_State{D} <: Abstract_State
    dof_indices::NTuple{D,Int} # tuple of indices to each physical dof. The order of each dof should be managed by the user
end
"add a show method for `Single_Particle_State`"
Base.show(io::IO, state::T) where {T<:Single_Particle_State} = print(io, "|$(state.dof_indices)⟩")


# add comparison between `Single_Particle_State`
Base.:(==)(state1::T, state2::T) where {T<:Single_Particle_State} = state1.dof_indices == state2.dof_indices
Base.:>(state1::T, state2::T) where {T<:Single_Particle_State} = state1.dof_indices > state2.dof_indices
Base.:<(state1::T, state2::T) where {T<:Single_Particle_State} = state1.dof_indices < state2.dof_indices

Base.isequal(state1::T, state2::T) where {T<:Single_Particle_State} = state1 == state2 # for sorting
Base.isless(state1::T, state2::T) where {T<:Single_Particle_State} = state1 < state2 # for sorting
Base.isgreater(state1::T, state2::T) where {T<:Single_Particle_State} = state1 > state2 # for sorting



"""
Highly-optimized In-place _Bubble Sort_ of the Many-particle State as a `Vector{<:Single_Particle_State}` with Returning of the Signature of Either `0` or a Boolean Value
---
Here `0` means the many-particle state is unphysical due to double occupation of single-particle states (if `is_hard_core=true`). While the boolean value means the parity of the legal many-particle state: `false` for `-1` and `true` for `+1`. This signature is crucial in constructing the second-quantized scattering amplitude, or the ED Hamiltonian.
- Args:
    - `ψ_many::Vector{<:Single_Particle_State}`: vector of single-particle states
    - `n_filled::Int`: number of filled single-particle states within the many-particle state, equal to `length(ψ_many)` (we explicitly pass the length to avoid the overhead of multiple call of `length()`)
- Named Args:
    - `quantum_statistics::Quantum_Statistics`: quantum statistics of the many-particle state, either `Bosonic` or `Fermionic`
    - `is_hard_core::Bool=true`: whether to ignore the unphysical many-body state due to double occupation when `is_hard_core=true`
---
For a small number of single-particle states (less than 30, for example), naive bubble sort, or insertion sort are much faster than any version of quick sort algorithm (due to zero overhead and zero memory allocation). 

Here it is tempting to use an integer `signature = +1,0,-1` to represent the parity of the many-particle state, but it turns out the even flip the sign of the integer `signature` can be much slower (cost about 220% more time) than simple bit flip of a bool value. So instead, we use a `Bool` value to represent the parity of the many-particle state, while still return 0 if there are duplicate single-particle states. Importantly, macro `@fastmath` can further speed up the bit flip operation by about 50% here! Slight improvement also comes from the dressing macro `@simd` and `@inbounds`.
"""
function sort_inplace_with_sign_flipped_fast!(ψ_many::Vector{S}, n_filled::Int; quantum_statistics::Quantum_Statistics, is_hard_core::Bool=true)::Union{Int,Bool} where {S<:Single_Particle_State}
    is_sign_flipped = false
    @simd for istate in 1:(n_filled-1) # `@simd` macro can improve about 20% performance
        @simd for iswap in 1:(n_filled-istate)
            @inbounds s1 = ψ_many[iswap]
            @inbounds s2 = ψ_many[iswap+1]
            if s1 == s2
                @match is_hard_core begin
                    true => begin
                        return 0 # unphysical state due to double occupation
                    end
                    false => nothing
                end
            elseif s1 > s2
                @inbounds ψ_many[iswap] = s2
                @inbounds ψ_many[iswap+1] = s1

                @match quantum_statistics begin
                    ::HilbertSpace.Bosonic => nothing
                    ::HilbertSpace.Fermionic => begin
                        @fastmath is_sign_flipped = ~is_sign_flipped # bit flip with macro `@fastmath` here can be 50% faster than using an integer signature and flip with `signature = -signature`
                    end
                end
            end
        end
    end
    return is_sign_flipped
end



"add a method without telling the length of the many-particle state as a `Vector{<:Single_Particle_State}`. This can be used for general situation"
function sort_inplace_with_sign_flipped_fast!(ψ_many::Vector{<:Single_Particle_State}; quantum_statistics::Quantum_Statistics, is_hard_core::Bool=true)::Union{Int,Bool}
    n_filled = length(ψ_many)
    if n_filled > 64
        @warn "The number of single-particle states is too large, probably quick sort algorithm can be better than bubble sort algorithm implemented for `sort_inplace_with_sign_flipped_fast!()`"
    end
    return sort_inplace_with_sign_flipped_fast!(ψ_many, n_filled; quantum_statistics=quantum_statistics, is_hard_core=is_hard_core)
end


"""
Highly-optimized In-place _Bubble Sort_ of the Many-particle State as a `Vector{<:Single_Particle_State}` _Without_ Returning the Signature 
---
The implementation here is just a simple wrapper of `sort_inplace_with_sign_flipped_fast!()`, but without returning the signature. This is useful when we only need to sort the many-particle state, but do not need to know the signature.
- Args:
    - `ψ_many::Vector{<:Single_Particle_State}`: vector of single-particle states
- Named Args:
    - `quantum_statistics::Quantum_Statistics`: quantum statistics of the many-particle state, either `Bosonic` or `Fermionic`
    - `is_hard_core::Bool=true`: whether to ignore the unphysical many-body state due to double occupation when `is_hard_core=true`
"""
function sort_inplace_fast!(ψ_many::Vector{<:Single_Particle_State}; quantum_statistics::Quantum_Statistics, is_hard_core::Bool=true)
    n_filled = length(ψ_many)
    signature = sort_inplace_with_sign_flipped_fast!(ψ_many, n_filled; quantum_statistics=quantum_statistics, is_hard_core=is_hard_core)

    return nothing
end



# TODO: add arithmic operations for `Single_Particle_State` to construct (symbolic) polynomials