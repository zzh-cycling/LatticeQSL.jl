function rdm_PXP(::Type{T}, subsystems::Vector{Int64}, state::Vector{ET}, pbc::Bool=true) where {N,T <: BitStr{N}, ET}
    # Usually subsystem indices count from the right of binary string. But this version we can count from the left, which is consistent with our intuition. This means that the systems we want to keep.
    # The function is to take common environment parts of the total basis, get the index of system parts in reduced basis, and then calculate the reduced density matrix.
    unsorted_basis = PXP_basis(T, pbc)
    @assert length(unsorted_basis) == length(state) "state length is expected to be $(length(unsorted_basis)), but got $(length(state))"
    
    subsystems=connected_components(subsystems)
    lengthlis=length.(subsystems)
    subsystems=vcat(subsystems...)
    mask = bmask(T, (N .-subsystems .+1)...)

    
    order = sortperm(unsorted_basis, by = x -> (takeenviron(x, mask), takesystem(x, mask))) #first sort by environment, then by system. The order of environment doesn't matter.
    basis, state = unsorted_basis[order], state[order]
    
    reduced_basis = move_subsystem.(T, joint_pxp_basis(lengthlis), Ref(subsystems))
    len = length(reduced_basis)
    # Initialize the reduced density matrix
    reduced_dm = zeros(ET, (len, len))

    # Keep track of indices where the key changes
    result_indices = Int[]
    current_key = -1
    for (idx, i) in enumerate(basis)
        key = takeenviron(i, mask)  # Get environment l bits
        if key != current_key
            @assert key > current_key "key is expected to be greater than $current_key, but got $key"
            push!(result_indices, idx)
            current_key = key
        end
    end
    # Add the final index to get complete ranges
    push!(result_indices, length(basis) + 1)

    for i in 1:length(result_indices)-1
        range = result_indices[i]:result_indices[i+1]-1         
        # Get indices in the reduced basis
        indices = searchsortedfirst.(Ref(reduced_basis), takesystem.(basis[range], mask))
        view(reduced_dm, indices, indices) .+= view(state, range) .* view(state, range)'
    end

    return reduced_dm
end
rdm_PXP(N::Int, subsystems::Vector{Int64}, state::Vector{ET}, pbc::Bool=true) where {ET} = rdm_PXP(BitStr{N, Int}, subsystems, state, pbc)


"""
    myprint(io::IO, xs...)

Enhanced print function with automatic flushing.

Prints arguments to the specified IO stream with an extra newline and
flushes the buffer to ensure immediate output, useful for real-time monitoring.

# Arguments
- `io::IO`: Output stream
- `xs...`: Arguments to print

# Example
```julia
myprint(stdout, "Results:", value1, value2)
```
"""
function myprint(io::IO, xs...)
    println(io, xs..., '\n')
    flush(io)
end

function ee(subrm::Matrix{ET}) where {ET}
    #  subrm=qi.ptrace(state*state',[2 for i in 1:N],[i for i in l+1:N])
    @assert ishermitian(subrm) "The reduced density matrix is not hermitian."
    spectrum=eigvals(subrm)
    EE=0
    for i in eachindex(spectrum)
        v=abs(spectrum[i])
            if v>1e-8
                EE+=-v*log(v)
            end
    end

    return EE
end