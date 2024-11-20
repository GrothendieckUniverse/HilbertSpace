#=====================================================================================================================================#
#======================================================== Single-Particle State ======================================================#
#=====================================================================================================================================#
"""
General Single-Particle State `Single_Particle_State{D}` <: Abstract_State`
---
with `D<:Integer` the number of physical d.o.f.
- Fields:
    - `dof_indices::Tuple`: tuple of indices to each physical d.o.f
"""
struct Single_Particle_State{D} <: Abstract_State where {D<:Integer}
    dof_indices::NTuple{D,Integer} # tuple of indices to each physical d.o.f
end
"add a show method for `Single_Particle_State`"
Base.show(io::IO, state::Single_Particle_State) = print(io, "ϕ_$(state.dof_indices)")


# add comparison between `Single_Particle_State`
Base.:(==)(state1::Single_Particle_State, state2::Single_Particle_State) = state1.dof_indices == state2.dof_indices
Base.:>(state1::Single_Particle_State, state2::Single_Particle_State) = state1.dof_indices > state2.dof_indices
Base.:<(state1::Single_Particle_State, state2::Single_Particle_State) = state1.dof_indices < state2.dof_indices


# TODO: add arithmic operations for `Single_Particle_State` to construct (symbolic) polynomials




#=====================================================================================================================================#
#========================================================= Many-Particle State =======================================================#
#=====================================================================================================================================#
"""
General Many Particle State `Many_Particle_State{N} <: Abstract_State`
---
with `N<:Integer` the number of single-particle states.
- Fields:
    - `single_particle_states::SVector{N,<:Single_Particle_State}`: static vector of single-particle states
"""
struct Many_Particle_State{N} <: Abstract_State where {N<:Integer}
    single_particle_states::SVector{N,<:Single_Particle_State}
end
"add a show method for `Many_Particle_State{N}`"
Base.show(io::IO, state::Many_Particle_State{N}) where {N} =
    for (i, ψ_single) in enumerate(state.single_particle_states)
        print(io, "$ψ_single")
        if i != N
            print(io, " ⊗ ")
        end
    end

"constructor of `Many_Particle_State` with a vector of single-particle states"
Many_Particle_State(single_particle_states::Vector{<:Single_Particle_State}) = Many_Particle_State(SVector{length(single_particle_states)}(single_particle_states))
"constructor of `Many_Particle_State` with a tuple of single-particle states"
Many_Particle_State(single_particle_states::NTuple{N,<:Single_Particle_State}) where {N} = Many_Particle_State(collect(single_particle_states))