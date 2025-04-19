abstract type Abstract_Hilbert_Space end


struct Single_Particle_Hilbert_Space <: Abstract_Hilbert_Space
    # ndof::Int # number of dof for the single-particle state
    # dof_name_to_range_pair_list::Vector{<:Pair{String,<:Vector{Int}}} # list of pair `dof_name::String => dof_range::Vector{Int}`. Here we want to keep the input order of each pair and the constructed single-particle state. So instead of using dictionary, here we choose to use vector of pairs
    # dof_name_list::Vector{String} # list of dof names (the order of the dof names is consistent with each single-particle state)

    state_list::Vector{<:Abstract_State} # list of single-particle states
    nstate::Int # number of single-particle states
    state_to_index_map::Dict{<:Abstract_State,Int} # hashmap `single-particle state -> index`

    metadata::Union{NamedTuple,Nothing} # metadata for the single-particle Hilbert space
end

"""
Constructor of `Single_Particle_Hilbert_Space`
---
- Args: 
    - `state_list::Vector{<:Abstract_State}`: list of single-particle states
- Named Args:
    - `metadata::Union{NamedTuple,Nothing}=nothing`: metadata for the single-particle Hilbert space
"""
function Single_Particle_Hilbert_Space(
    state_list::Vector{<:Abstract_State};
    metadata::Union{NamedTuple,Nothing}=nothing
)
    nstate = length(state_list)
    state_to_index_map = Dict(
        state => i
        for (i, state) in enumerate(state_list)
    )

    return Single_Particle_Hilbert_Space(
        state_list,
        nstate,
        state_to_index_map,
        metadata,
    )
end


"""
Another Constructor of `Single_Particle_Hilbert_Space`
---
with `dof_name_to_range_pair_list` and `metadata`.
- Named Args:
    - `dof_name_to_range_pair_list::Dict{String,<:Vector{Int}}`: list of pair `dof_name::String => dof_range::Vector{Int}`
"""
function Single_Particle_Hilbert_Space(;
    dof_name_to_range_pair_list::Vector{<:Pair{String,<:Vector{Int}}},
    metadata::Union{NamedTuple,Nothing}=nothing,
)
    ndof = length(dof_name_to_range_pair_list)
    dof_range_list = last.(dof_name_to_range_pair_list)
    dof_name_list = first.(dof_name_to_range_pair_list)

    single_particle_state_list = [Single_Particle_State(dof_tuple) for dof_tuple in Iterators.product(dof_range_list...)] |> vec # force it to be vector of single-particle states

    nstate = length(single_particle_state_list)
    single_particle_state_to_index_map = Dict(
        state => i
        for (i, state) in enumerate(single_particle_state_list)
    )



    return Single_Particle_Hilbert_Space(
        # ndof,
        # dof_name_to_range_pair_list,
        # dof_name_list,
        single_particle_state_list,
        nstate,
        single_particle_state_to_index_map,
        metadata
    )
end



# TODO: `Infinite_Single_Particle_Hilbert_Space` for infinite-dimensional single-particle Hilbert space



"""
Struct `Many_Particle_Hilbert_Space <: `
---
- Fields:
    - `state_list::Vector{<:Vector{<:Abstract_State}}`: Here we explicitly store the entire list of the many-particle state for fastest iteration performance (required for ED Hamiltonian construction, for example), although it may cost HUGE memory allocations here!. Note: in principle one can avoid such huge memory allocation with the lazy generator. However, the cost is that the iteration speed will be more than 100 times slower (I have test this!) due to overhead of method call, internal memory allocation within each step (for complicated generator), as well as the cache missing (for simple generator)
    - `nstate::Int`: number of many-body states
    - `state_to_index_map::Dict{<:Vector{<:Abstract_State},Int}`: hashmap `many-particle state -> index`
"""
struct Many_Particle_Hilbert_Space <: Abstract_Hilbert_Space
    nparticle::Int # number of particles
    state_list::Vector{<:Vector{<:Abstract_State}} # Here we explicitly store the entire list of the many-particle state for fastest iteration performance (required for ED Hamiltonian construction, for example), although it may cost HUGE memory allocations here!. Note: in principle one can avoid such huge memory allocation with the lazy generator. However, the cost is that the iteration speed will be more than 100 times slower (I have test this!) due to overhead of method call, internal memory allocation within each step (for complicated generator), as well as the cache missing (for simple generator)
    nstate::Int # number of many-body states (which can be a large number so it's better to store it explicitly)
    state_to_index_map::Dict{<:Vector{<:Abstract_State},Int} # hashmap `many-particle state -> index`
end


"""
Constructor of `Many_Particle_Hilbert_Space`
---
- Args:
    - `many_particle_state_list::Vector{<:Vector{<:Abstract_State}}`: list of many-particle states
"""
function Many_Particle_Hilbert_Space(;
    many_particle_state_list::Vector{<:Vector{<:Abstract_State}}
)
    nparticle = length(many_particle_state_list[1]) # number of particles
    # check if all many-particle states have the same number of particles
    @assert all(length(many_particle_state) == nparticle for many_particle_state in many_particle_state_list) "All many-particle states should have the same number of particles"

    nstate = length(many_particle_state_list)

    many_particle_state_to_index_map = Dict(
        many_particle_state => i
        for (i, many_particle_state) in enumerate(many_particle_state_list)
    )

    return Many_Particle_Hilbert_Space(
        nparticle,
        many_particle_state_list,
        nstate,
        many_particle_state_to_index_map
    )
end




# TODO: `Infinite_Many_Particle_Hilbert_Space` for infinite-dimensional many-particle Hilbert space
