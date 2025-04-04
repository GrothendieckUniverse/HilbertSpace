
abstract type Finite_Dimensional_Hilbert_Space <: Abstract_Hilbert_Space end


#=====================================================================================================================================#
#=================================================== Single-Particle Hilbert Space ===================================================#
#=====================================================================================================================================#
"""
Struct `Finite_Dimensional_Single_Particle_Hilbert_Space <: Finite_Dimensional_Hilbert_Space`
---
- Fields:
    - `ndof`: number of d.o.f for the single-particle state
    - `dof_name`: name of each d.o.f
    - `dof_range_map`: hashmap from dof name to the index range of each d.o.f
    - `state_list`: list of single-particle states
    - `nstate`: number of single-particle states
    - `state_to_index_map`: hashmap `single-particle state -> index`
"""
struct Finite_Dimensional_Single_Particle_Hilbert_Space <: Finite_Dimensional_Hilbert_Space
    ndof::Int # number of d.o.f for the single-particle state
    dof_range_map::Dict{String,UnitRange{Int}} # hashmap from dof name to the index range of each d.o.f

    state_list::Vector{<:Single_Particle_State} # list of single-particle states
    nstate::Int # number of single-particle states
    state_to_index_map::Dict{<:Single_Particle_State,Int} # hashmap `single-particle state -> index`
end

"""
Constructor of `Finite_Dimensional_Single_Particle_Hilbert_Space`
---
- Named Args:
    - `dof_range_map::Dict{String,UnitRange{Int}}`: hashmap from dof name to the index range of each d.o.f
"""
function Finite_Dimensional_Single_Particle_Hilbert_Space(;
    dof_range_map::Dict{String,UnitRange{Int}},
)
    ndof = length(keys(dof_range_map))
    dof_range_list = values(dof_range_map) |> collect
    dof_name_list = keys(dof_range_map) |> collect

    single_particle_state_list = [Single_Particle_State(finite_dof_indices) for finite_dof_indices in Iterators.product(dof_range_list...)] |> vec # force it to be vector of single-particle states

    nstate = length(single_particle_state_list)
    single_particle_state_to_index_map = Dict(
        state => i
        for (i, state) in enumerate(single_particle_state_list)
    )

    return Finite_Dimensional_Single_Particle_Hilbert_Space(
        ndof,
        dof_range_map,
        single_particle_state_list,
        nstate,
        single_particle_state_to_index_map,
    )
end



# """
# Struct `Finite_Dimensional_Multi_Particle_Hilbert_Space <: Finite_Dimensional_Hilbert_Space`
# ---
# - Fields:
#     - `ndof`: number of d.o.f for the single-particle state
#     - `dof_range_map`: hashmap from dof name to the index range of each d.o.f for the single-particle state
#     - `nparticle`: number of particles
#     - `state_list`: list of multi-particle states. Here we would like to arrange each multi-particle state `(s1,s2,...)` in a way such that `s1<=s2<=...` to save memory (this may also be helpful for fermionic parity if necessary)
#     - `nstate`: number of multi-particle states
#     - `state_to_index_map`: hashmap `multi-particle state -> index`
# """
# struct Finite_Dimensional_Multi_Particle_Hilbert_Space{D} <: Finite_Dimensional_Hilbert_Space where {D<:Int}
#     nparticle::Int # number of particles
#     state_list::Vector{<:NTuple{D,<:Single_Particle_State}} # list of multi-particle states. Here we would like to arrange each multi-particle state `(s1,s2,...)` in a way such that `s1<=s2<=...` to save memory (this may also be helpful for fermionic parity if necessary)
#     nstate::Int # number of multi-particle states
#     state_to_index_map::Dict{<:NTuple{D,<:Single_Particle_State},<:Int} # hashmap `multi-particle state -> index`
# end



# TODO: `Infinite_Single_Particle_Hilbert_Space` for infinite-dimensional single-particle Hilbert space





#=====================================================================================================================================#
#==================================================== Many-Particle Hilbert Space ====================================================#
#=====================================================================================================================================#

"""
Struct `Finite_Dimensional_Many_Particle_Hilbert_Space <: Finite_Dimensional_Hilbert_Space`
---
- Fields:
    - `state_list::Vector{<:Vector{<:Single_Particle_State}}`: Here we explicitly store the entire list of the many-particle state for fastest iteration performance (required for ED Hamiltonian construction, for example), although it may cost HUGE memory allocations here!. Note: in principle one can avoid such huge memory allocation with the lazy generator. However, the cost is that the iteration speed will be more than 100 times slower (I have test this!) due to overhead of method call, internal memory allocation within each step (for complicated generator), as well as the cache missing (for simple generator)
    - `nstate::Int`: number of many-body states
    - `state_to_index_map::Dict{<:Vector{<:Single_Particle_State},Int}`: hashmap `many-particle state -> index`
"""
struct Finite_Dimensional_Many_Particle_Hilbert_Space <: Finite_Dimensional_Hilbert_Space
    nparticle::Int # number of particles
    state_list::Vector{<:Vector{<:Single_Particle_State}} # Here we explicitly store the entire list of the many-particle state for fastest iteration performance (required for ED Hamiltonian construction, for example), although it may cost HUGE memory allocations here!. Note: in principle one can avoid such huge memory allocation with the lazy generator. However, the cost is that the iteration speed will be more than 100 times slower (I have test this!) due to overhead of method call, internal memory allocation within each step (for complicated generator), as well as the cache missing (for simple generator)
    nstate::Int # number of many-body states (which can be a large number so it's better to store it explicitly)
    state_to_index_map::Dict{<:Vector{<:Single_Particle_State},Int} # hashmap `many-particle state -> index`
end


"""
Constructor of `Finite_Dimensional_Many_Particle_Hilbert_Space`
---
- Args:
    - `many_particle_state_list::Vector{<:Vector{<:Single_Particle_State}}`: list of many-particle states
"""
function Finite_Dimensional_Many_Particle_Hilbert_Space(;
    many_particle_state_list::Vector{<:Vector{<:Single_Particle_State}}
)
    nparticle = length(many_particle_state_list[1]) # number of particles
    # check if all many-particle states have the same number of particles
    @assert all(length(many_particle_state) == nparticle for many_particle_state in many_particle_state_list) "All many-particle states should have the same number of particles"

    nstate = length(many_particle_state_list)

    many_particle_state_to_index_map = Dict(
        many_particle_state => i
        for (i, many_particle_state) in enumerate(many_particle_state_list)
    )

    return Finite_Dimensional_Many_Particle_Hilbert_Space(
        nparticle,
        many_particle_state_list,
        nstate,
        many_particle_state_to_index_map
    )
end




# TODO: `Infinite_Many_Particle_Hilbert_Space` for infinite-dimensional many-particle Hilbert space
