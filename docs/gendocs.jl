#!/usr/bin/env julia

module_name = :TTCal
output_external = "external.md"
output_internal = "internal.md"

@eval using $module_name
import Base.Markdown.plain
import Base.Docs: FuncDoc, TypeDoc

cd(dirname(@__FILE__))

exports = @eval names($module_name)
meta    = @eval Docs.meta($module_name)

# Get the plain text docs corresponding to each type and function.
docs = Dict{Any,UTF8String}()

function add_docs(obj,doc::FuncDoc)
    md = collect(values(doc.meta))
    docs[symbol(obj)] = plain(md)
end
function add_docs(obj,doc::TypeDoc)
    md  = doc.main
    str = plain(md)
    # search for all the constructors
    for signature in doc.order
        str = string(str,"\n",plain(doc.meta[signature]))
    end
    docs[symbol(obj)] = str
end
add_docs(obj,doc) = nothing

for (obj,doc) in meta
    add_docs(obj,doc)
end

# Check to see if any exported symbols are missing documentation.
for sym in exports
    sym == module_name && continue
    sym = symbol(module_name,".",sym)
    if !haskey(docs,sym)
        warn("$sym is exported and does not have any documentation.")
        docs[sym] = "\n"
    end
end

# Collect the unexported symbols with documentation
not_exports = Symbol[]
for sym in keys(docs)
    str = replace(string(sym),"$module_name.","")
    stripped_sym = symbol(str)
    stripped_sym in exports && continue
    push!(not_exports,symbol(str))
end
sort!(not_exports)

# Write the documentation
open(output_external,"w") do file
    write(file,"""
        <!---
        This is an auto-generated file and should not be edited directly.
        -->

        ## API

        """)
    for sym in exports
        sym == module_name && continue
        write(file,"### $sym\n\n")
        sym = symbol(module_name,".",sym)
        write(file,docs[sym])
        write(file,"\n")
    end
end
open(output_internal,"w") do file
    write(file,"""
        <!---
        This is an auto-generated file and should not be edited directly.
        -->

        ## Internal Documentation

        """)
    for sym in not_exports
        write(file,"### $sym\n\n")
        sym = symbol(module_name,".",sym)
        write(file,docs[sym])
        write(file,"\n")
    end
end

