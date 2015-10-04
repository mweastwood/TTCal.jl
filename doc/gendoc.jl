#!/usr/bin/env julia

module_name = :TTCal
output = "api.md"

@eval using $module_name
import Base.Markdown.plain
import Base.Docs: FuncDoc, TypeDoc

cd(dirname(@__FILE__))

exports = names(TTCal)
meta    = Docs.meta(TTCal)

# Get the plain text docs corresponding to each type and function.
docs = Dict{Any,UTF8String}()

function add_docs(obj,doc::FuncDoc)
    md = collect(values(doc.meta))
    docs[symbol(obj)] = plain(md)
end
function add_docs(obj,doc::TypeDoc)
    md = doc.main
    docs[symbol(obj)] = plain(md)
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

# Write the documentation
open(output,"w") do file
    write(file,"""
        <!---
        This is an auto-generated file and should not be edited directly.
        -->

        """)
    for sym in exports
        sym == module_name && continue
        write(file,"## $sym\n\n")
        sym = symbol(module_name,".",sym)
        write(file,docs[sym])
    end
end

