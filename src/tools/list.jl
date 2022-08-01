#import Base.show, Base.push!, Base.delete!, Base.length

abstract type AbstractNode end

mutable struct Node{T}
    data::T
    next::Node{T}
    prev::Node{T}
    function Node{T}() where T
        return new()
    end
    function Node{T}(data::T) where T
        return new(data)
    end
end

function Base.show(io::IO, node::Node{T}) where T
    print(io, "Node(data: $(node.data))")
end

function Base.getindex(node::Node, index::Int)
    nodei = node
    if index>0
        for i in 1:index
            nodei = nodei.next
        end
    end

    if index<0
        index = abs(index)
        for i in 1:index
            nodei = nodei.prev
        end
    end

    return nodei
end

#function Node{T}(data::T)
    #@show 1
    #node = Node{T}()
    #@show 1
    #node.data = data
#end

export List
mutable struct List{T}
    size  ::Int
    first::Node{T}
    last ::Node{T}
    function List{T}() where T
        return new(0)
    end
    function List{T}(arr::Array{T,1}) where T
        list = new(0)
        length(arr)==0 && return

        for obj in arr
            node = Node{T}(obj)
            if list.size==0
                list.first= node
                list.last = node
                node.prev = node
                node.next = node
            else
                node.prev = list.last
                node.next = list.first # circular
                list.first.prev = node
                list.last.next  = node
                list.last = node
            end
            list.size += 1
        end
        return list
    end
end

function Base.length(list::List)
    return list.size
end


#type Iterator{T}
    #container::T
#end

#iterator(container) = Iterator(container)

# DummyNode used as start point for iteration
#type DummyNode{T}
    #next::Node{T}
#end

#Base.start{T}(list::List{T}) = DummyNode{T}(list.first)
#Base.next{T}(list::List{T}, node) = node.next, node.next
#Base.done{T}(list::List{T}, node) = node == list.last

function Base.iterate(list::List{T}, state=(nothing,0)) where T
    node, count = state
    if node===nothing
        return (list.first, (list.first.next, 1))
    elseif node.prev!=list.last
        return (node, (node.next, count+1))
    else
        return nothing
    end
end

export check
function check(list::List{T}) where T
    node = list.first
    for i in 1:length(list)
        if objectid(node.next.prev)==objectid(node)
            println("ok")
        else
            println("error")
        end
        node = node.next
    end
end


function Base.push!(list::List{T}, node::Node{T}) where T
    if list.size==0
        list.first= node
        list.last = node
        node.prev = node
        node.next = node
    else
        node.prev = list.last
        node.next = list.first # circular
        list.first.prev = node
        list.last.next  = node
        list.last = node
    end
    list.size += 1
end

function Base.push!(list::List{T}, data::T) where T
    push!(list, Node{T}(data))
end

function Base.insert!(list::List{T}, nodepos::Node{T}, newnode::Node{T}) where T
    @assert list.size>0

    list.first==nodepos && (list.first = newnode)

    newnode.prev = nodepos.prev
    newnode.prev.next = newnode
    newnode.next = nodepos
    nodepos.prev = newnode

    list.size += 1
end

function Base.delete!(list::List{T}, node::Node{T}) where T
    @assert list.size > 0

    list.first==node && (list.first = node.next)
    list.last ==node && (list.last  = node.prev)

    node.prev.next = node.next
    node.next.prev = node.prev
    list.size -= 1
end

#=

l = List{Int}()
n = Node{Int}(1)
push!(l, n)
n = Node{Int}(2)
push!(l, n)
n = Node{Int}(3)
push!(l, n)
n = Node{Int}(4)
push!(l, n)
ii = n
n = Node{Int}(44)
insert!(l, ii, n)

delete!(l, n)

for node in l
    @show "========="
    @show node.prev.data
    @show node.data
    @show node.next.data
    #@show (node.data, node.prev, node.next)
end

node = l.first
@show node[1]
@show node[2]
@show node[3]
@show node[4]
node = l.last
@show node[-1]
@show node[-2]
@show node[-3]

#A*B = Aik*Bkj
#A*B*C = Aik*Bkm*Cmj

=#
