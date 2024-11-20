
abstract type Finite_Dimensional_Hilbert_Space <: Abstract_Hilbert_Space end


#=====================================================================================================================================#
#=================================================== Single-Particle Hilbert Space ===================================================#
#=====================================================================================================================================#
"""
Struct `Finite_Dimensional_Single_Particle_Hilbert_Space{D} <: Finite_Dimensional_Hilbert_Space`
---
with `D<:Integer` the number of indices to label the single-particle states.
- Fields:
    - `ndof::Integer`: number of d.o.f for the single-particle state: always explicitly stored as the parametric `D<:Integer` (please distinguish from `dof_ndof` below!)
    - `dof_name::SVector{D,String}`: static vector of the names for each d.o.f of the single-particle state
    - `dof_ndof::SVector{D,<:Integer}`: static vector of the number of states for each d.o.f of the single-particle state. For example, `dof_ndof = (2, 3)` means the first d.o.f has 2 states and the second has 3 states
    - `state_list::Vector{<:Single_Particle_State}`: list of single-particle states
    - `nstate::Int`: number of single-particle states
    - `state_to_index_map::Dict{<:Single_Particle_State,<:Integer}`: hashmap `single-particle state -> index`
"""
struct Finite_Dimensional_Single_Particle_Hilbert_Space{D} <: Finite_Dimensional_Hilbert_Space where {D<:Integer}
    ndof::Integer # number of d.o.f for the single-particle state: always explicitly stored as the parametric `D<:Integer` (please distinguish from `dof_ndof` below!)
    dof_name::SVector{D,String} # static vector of the names for each d.o.f of the single-particle state
    dof_ndof::SVector{D,<:Integer} # static vector of the number of states for each d.o.f of the single-particle state. For example, `dof_ndof = (2, 3)` means the first d.o.f has 2 states and the second has 3 states

    state_list::Vector{<:Single_Particle_State} # list of single-particle states
    nstate::Int # number of single-particle states
    state_to_index_map::Dict{<:Single_Particle_State,<:Integer} # hashmap `single-particle state -> index`
end

"""
Constructor of `Finite_Dimensional_Single_Particle_Hilbert_Space`
---
- Named Args:
    - `dof_ndof::NTuple{D,<:Integer}`: vector of the number of degrees of freedom for each d.o.f of the single-particle state. For example, `dof_ndof = (2, 3)` means the first d.o.f has 2 states and the second has 3 states
    - `dof_name::NTuple{D,String}`: vector of the name for each d.o.f of the single-particle state
"""
function Finite_Dimensional_Single_Particle_Hilbert_Space(;
    dof_ndof::Vector{<:Integer},
    dof_name::Vector{String}=Vector{String}(),
)
    D = length(dof_ndof)
    if isempty(dof_name)
        dof_name = SVector{D}(["dof_$i" for i in 1:D])
    end

    single_particle_state_list = [Single_Particle_State(finite_dof_indices) for finite_dof_indices in Iterators.product((1:n for n in dof_ndof)...)] |> vec # force it to be vector of single-particle states

    nstate = length(single_particle_state_list)
    single_particle_state_to_index_map = Dict(
        state => i
        for (i, state) in enumerate(single_particle_state_list)
    )

    return Finite_Dimensional_Single_Particle_Hilbert_Space{D}(
        D,
        SVector{D}(dof_name),
        SVector{D}(dof_ndof),
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
Struct `Finite_Dimensional_Many_Particle_Hilbert_Space{N,I} <: Finite_Dimensional_Hilbert_Space`
---
with `N<:Integer` the number of single-particle states for each many-particle state, i.e. the number of particles of the many-particle system, and `I<:Union{Int,BigInt}` the type of the number of many-body states.
- Fields:
    - `state_generator::Base.Generator{<:Combinatorics.Combinations,<:Function}`: lazy generator for the list of many-body states. Here each many-body state is a vector of single-particle states. Note: we should never try to explicitly store the Full list of many-body states, since it can cost formidablely large amount of memory! On the other hand, it is suitable for the on-the-fly algorithm for ED calculation
    - `nstate::I`: number of many-body states (it can be a large so it's better to store it explicitly)
    - `state_to_index_map::Dict{<:Many_Particle_State{N},I}`: hashmap `many-particle state -> index`
"""
struct Finite_Dimensional_Many_Particle_Hilbert_Space{N,I} <: Finite_Dimensional_Hilbert_Space where {N<:Integer,I<:Union{Int,BigInt}}
    state_generator::Base.Generator{<:Combinatorics.Combinations,<:Function} # lazy generator for the list of many-body states. Here each many-body state is a vector of single-particle states. Note: we should never try to explicitly store the Full list of many-body states, since it can cost formidablely large amount of memory! On the other hand, it is suitable for the on-the-fly algorithm for ED calculation
    nstate::I # number of many-body states (it can be a large so it's better to store it explicitly)
    state_to_index_map::Dict{<:Many_Particle_State{N},I} # hashmap `many-particle state -> index`
end

"""
Constructor of `Finite_Dimensional_Many_Particle_Hilbert_Space`
---
- Named Args:
    - `single_particle_Hilbert_space::Finite_Dimensional_Single_Particle_Hilbert_Space`: single-particle Hilbert space
    - `nfermion::N`: number of fermions
"""
function Finite_Dimensional_Many_Particle_Hilbert_Space(;
    single_particle_Hilbert_space::Finite_Dimensional_Single_Particle_Hilbert_Space,
    nfermion::Int
)
    g = Combinatorics.combinations(single_particle_Hilbert_space.state_list, nfermion)

    many_particle_state_generator = (g.f.reorder.a[v] for v in g.iter)

    nstate = let nstate_big::BigInt = binomial(BigInt(single_particle_Hilbert_space.nstate), nfermion) # convert to BigInt to avoid overflow
        nstate_big <= typemax(Int) ? Int(nstate_big) : nstate_big # shrink to Integer if possible
    end
    nstate_type = typeof(nstate)

    many_particle_state_to_index_map = Dict(
        Many_Particle_State(many_particle_state) => nstate_type(i)
        for (i, many_particle_state) in enumerate(many_particle_state_generator)
    )

    return Finite_Dimensional_Many_Particle_Hilbert_Space{nfermion,nstate_type}(
        many_particle_state_generator,
        nstate,
        many_particle_state_to_index_map
    )
end


# TODO: `Infinite_Many_Particle_Hilbert_Space` for infinite-dimensional many-particle Hilbert space
