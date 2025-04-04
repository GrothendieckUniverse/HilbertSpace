
abstract type Finite_Dimensional_Hilbert_Space <: Abstract_Hilbert_Space end


#=====================================================================================================================================#
#=================================================== Single-Particle Hilbert Space ===================================================#
#=====================================================================================================================================#
"""
Struct `Finite_Dimensional_Single_Particle_Hilbert_Space <: Finite_Dimensional_Hilbert_Space`
---
- Fields:
    - `ndof`: number of dof for the single-particle state
    - `dof_name_to_range_pair_list`: list of pair `dof_name::String => dof_range::Vector{Int}`. Here we want to keep the input order of each pair and the constructed single-particle state. So instead of using dictionary, here we choose to use vector of pairs
    - `dof_name_list`: list of dof names (the order of the dof names is consistent with each single-particle state)
    - `state_list`: list of single-particle states
    - `nstate`: number of single-particle states
    - `state_to_index_map`: hashmap `single-particle state -> index`
"""
struct Finite_Dimensional_Single_Particle_Hilbert_Space <: Finite_Dimensional_Hilbert_Space
    ndof::Int # number of dof for the single-particle state
    dof_name_to_range_pair_list::Vector{<:Pair{String,<:Vector{Int}}} # list of pair `dof_name::String => dof_range::Vector{Int}`. Here we want to keep the input order of each pair and the constructed single-particle state. So instead of using dictionary, here we choose to use vector of pairs
    dof_name_list::Vector{String} # list of dof names (the order of the dof names is consistent with each single-particle state)
    state_list::Vector{<:Single_Particle_State} # list of single-particle states
    nstate::Int # number of single-particle states
    state_to_index_map::Dict{<:Single_Particle_State,Int} # hashmap `single-particle state -> index`
end

"""
Constructor of `Finite_Dimensional_Single_Particle_Hilbert_Space`
---
- Named Args:
    - `dof_name_to_range_pair_list::Dict{String,<:Vector{Int}}`: list of pair `dof_name::String => dof_range::Vector{Int}`
"""
function Finite_Dimensional_Single_Particle_Hilbert_Space(;
    dof_name_to_range_pair_list::Vector{<:Pair{String,<:Vector{Int}}},
)
    ndof = length(dof_name_to_range_pair_list)
    dof_range_list = last.(dof_name_to_range_pair_list)
    dof_name_list = first.(dof_name_to_range_pair_list)

    single_particle_state_list = [Single_Particle_State(finite_dof_indices) for finite_dof_indices in Iterators.product(dof_range_list...)] |> vec # force it to be vector of single-particle states

    nstate = length(single_particle_state_list)
    single_particle_state_to_index_map = Dict(
        state => i
        for (i, state) in enumerate(single_particle_state_list)
    )

    return Finite_Dimensional_Single_Particle_Hilbert_Space(
        ndof,
        dof_name_to_range_pair_list,
        dof_name_list,
        single_particle_state_list,
        nstate,
        single_particle_state_to_index_map,
    )
end



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
