
"""
findclosing(opening, closing, string, start)
"""
function findclosing(opening::AbstractString, closing::AbstractString, str::AbstractString, start::Int)
    @assert opening!=closing
    pattern = "(?:$opening|$closing)"
    if length(opening)==1 && length(closing)==1
        pattern = "(?<!\\)(?:$opening|$closing)"
    end
    #@show pattern
    pattern = replace(pattern, "\\" => "\\\\")
    #@show pattern
    #@show Regex(pattern)
    pos = start
    depth = 1
    while true
        rng = findnext(Regex(pattern), str, pos)
        rng==nothing && return nothing
        if str[rng]==opening
            depth += 1
        else
            depth -= 1
            depth==0 && return rng
        end
        pos = rng.stop+1
    end
end

function findclosure(opening::AbstractString, closing::AbstractString, str::AbstractString, start::Int)
    @assert opening!=closing
    openingpat = opening
    closingpat = closing

    if occursin(opening, "[]")
        openingpat = "\\"*opening
    end
    if occursin(closing, "[]")
        closingpat = "\\"*closing
    end

    pattern = "(?:$openingpat|$closingpat)"
    if length(opening)==1 && length(closing)==1
        pattern = "(?<!\\\\)(?:$openingpat|$closingpat)"
    end
    pos = start
    depth = 0
    while true

        rng = findnext(Regex(pattern), str, pos)
        rng==nothing && return nothing
        if str[rng]==opening
            depth==0 && (start=rng.start)
            depth += 1
        else
            depth -= 1
            depth==0 && return start:rng.stop
        end
        pos = rng.stop+1
    end
end

export tex2xml

function tex2xml(filename::String)
    #xdoc = Xdoc()
    root = Xnode("latex")
    text = read(filename, String)
    root = readtex(text, 1:latexindex(text))
    return root
end

function readtex(text::String, rng::UnitRange{Int})
    pardelim = r"^ \s* (?: \$\$   |
                    \\[][] |
                    \\ (?: b(?:egin   | ibitem)            |
                             c(?:aption | hapter)            |
                             end                             |
                             footnote                        |
                             item                            |
                             label                           |
                             marginpar                       |
                             n(?:ew(?:lin | pag)e | oindent) |
                             par(?:agraph | box | t)         |
                             s(?:ection |
                                 ub(?:paragraph |
                                      s(?:(?:ubs)?section)
                                   )
                              )                              |
                             [a-z]*(page | s(?:kip | space))
                          ) \b
                  )"mx

    children = Xnode[]

    len = length(text)
    pos = rng.start
    while i<=rng.stop
        pos = nextind(text, pos-1)
        c = text[pos]
        if !isprint(c) 
            continue
        end

        newline = c=='\n' || pos==1

        if newline
            if match(r"^\s*$"m, line)!=nothing # blank line
                if length(root.children)>0 && root.children[end].name!="blank"
                    xnode = Xnode("blank", Dict(), "")
                    push!(root.children, xnode)
                end
                pos = rng.stop
            elseif match(r"^\s*%.*?$"m, line)!=nothing # comment
                if length(root.children)>0 && root.children[end].name=="comment"
                    root.children[end].content *= "\n"*line
                else
                    xnode = Xnode("comment", Dict(), line)
                    push!(root.children, xnode)
                end
                pos = rng.stop
            end
            continue
        end

        if c=='\\' 
            if isletter(text[pos+1])
                m = match(r"(\\[a-zA-Z]*\b)"m, text, pos)
                pos += length(m.match)
                name = m.captures[1]

                if name=="begin"
                    m = match(r"{\\[a-zA-Z]*\b}"m, text, pos)
                    pos += length(m.match)
                    env = m.captures[1]
                    opening = "\\begin{$env}"
                    closing = "\\end{$env}"
                    rng = findclosure(opening, closing, text, pos)
                    xtex = readtex(text, rng.start+length(opening):rng.stop-length(closing))
                    xtex.name = "environment"
                    xtex.attributes["name"] = name
                    push!(children, xtex)
                    pos = rng.stop+1
                else
                    children = Xnode[]
                    while true
                        rng = findnext(r"\S", text, pos)
                        if text[rng] == "["
                            rng = findclosure("[", "]", text, rng.start)
                            rng==nothing && break
                            push!(children, Xnode("squarebraces", Dict(), text[rng]))
                            pos = rng.stop+1
                        elseif text[rng] == "{"
                            rng = findclosure("{", "}", text, rng.start)
                            rng==nothing && break
                            xtex = readtex(text, rng.start+1:rng.stop-1)
                            xtex.name = "curlybraces"
                            push!(children, xtex)
                            pos = rng.stop+1
                        #elseif text[rng] == "\$"
                            #rng = findnext(r"\$"m, text, pos)
                            #rng==nothing && error()
                            #xtex = readtex(text, pos:rng.stop-1)
                            #xnode = Xnode("inline-equation", Dict(), xtex.children)
                        else
                            break
                            #rng = findnext(r"(\[|{|\\)", text, pos)
                        end
                        #pos = rng.stop
                    end
                    xnode = Xnode("command", Dict("name"=>name))
                    push!(children, xnode)

                end

            end
            #rof = text[findnext(r".*?$", text, pos)]
            #m=match(r"(\\[a-zA-Z]*\b)"m, rof)
            #if m!=nothing && m.match.offset==0
            #end
        elseif c=='{' 
            rng = findclosure("{", "}", text, pos)
            rng==nothing && break
            xtex = readtex(text, rng.start+1:rng.stop-1)
            xnode = Xnode("curlybraces", Dict(), xtex.children)
            push!(children, xnode)
            pos = rng.stop+1

        else # raw text
            #rex = r"( \( | \) | { | } | \\ | %"x
            #rng = findnext(r"(\[|{|\\)", text, pos)
            rng = findnext(pardelim, text, pos)
            rng==nothing && break
            xtex = Xnode("paragraph", Dict(), [], text[pos:rng.start-1])
            pos = rng.start

        end


    end

    return Xnode("tex", Dict(), children)

end
function tex2xml2(filename::String)
    #xdoc = Xdoc()
    root = Xnode("latex")
    text = read(filename, String)

    # Preable
    # ≡≡≡≡≡≡≡≡≡
    
    pos = 1
    while true
        # current line
        rng = findnext(r"^.*?$"m, text, pos)
        rng==nothing && break
        line = text[rng]

        match(r"^\s*\\begin{document}", line)!=nothing && break

        if match(r"^\s*$"m, line)!=nothing # blank line
            if length(root.children)>0 && root.children[end].name!="blank"
                xnode = Xnode("blank", Dict(), "")
                push!(root.children, xnode)
            end
            pos = rng.stop+1
        elseif match(r"^\s*%.*?$"m, line)!=nothing # comment
            if length(root.children)>0 && root.children[end].name=="comment"
                root.children[end].content *= "\n"*line
            else
                xnode = Xnode("comment", Dict(), line)
                push!(root.children, xnode)
            end
            pos = rng.stop+1
        elseif (m=match(r"^\s*(\\[a-zA-Z]*\b)"m, line)) != nothing # command
            name = m.captures[1]
            xnode = Xnode("command", Dict("name"=>name))
            rng = findnext(r"^\s*\\[a-zA-Z]*\b"m, text, pos)
            pos = rng.stop+1

            while true
                rng = findnext(r"\S", text, pos)
                if text[rng] == "["
                    rng = findclosure("[", "]", text, rng.start)
                    rng==nothing && break
                    push!(xnode.children, Xnode("arg", Dict(), text[rng]))
                    pos = rng.stop+1
                elseif text[rng] == "{"
                    rng = findclosure("{", "}", text, rng.start)
                    rng==nothing && break
                    push!(xnode.children, Xnode("part", Dict(), text[rng]))
                    pos = rng.stop+1
                else
                    rng = findnext(r"$"m, text, pos)
                    pos = rng.stop+1
                    break
                end
                pos = rng.stop+1
            end

            push!(root.children, xnode)
        else
            error("something wrong")

        end

    end

    # Document
    
    linecommand = r"^ \s* (?: \$\$   |
                    \\[][] |
                    \\ (?: b(?:egin   | ibitem)            |
                             c(?:aption | hapter)            |
                             end                             |
                             footnote                        |
                             item                            |
                             label                           |
                             marginpar                       |
                             n(?:ew(?:lin | pag)e | oindent) |
                             par(?:agraph | box | t)         |
                             s(?:ection |
                                 ub(?:paragraph |
                                      s(?:(?:ubs)?section)
                                   )
                              )                              |
                             [a-z]*(page | s(?:kip | space))
                          ) \b
                  )"mx
    #newlinepat = replace(newlinepat, r"\s" => "")
    #newlineregex = Regex(newlinepat, "m")
    #@show newlineregex

    pos = 1
    while true
        # current line
        rng = findnext(r"^.*?$"m, text, pos)
        rng==nothing && break
        line = text[rng]
        #match(r"^\s*\\begin{end}", line)!=nothing && break

        if match(r"^\s*$"m, line)!=nothing # blank line
            if length(root.children)>0 && root.children[end].name!="blank"
                xnode = Xnode("blank", Dict(), "")
                push!(root.children, xnode)
            end
            pos = rng.stop+1
        elseif match(r"^\s*%.*?$"m, line)!=nothing # comment
            if length(root.children)>0 && root.children[end].name=="comment"
                root.children[end].content *= "\n"*line
            else
                xnode = Xnode("comment", Dict(), line)
                push!(root.children, xnode)
            end
            pos = rng.stop+1
        elseif (m=match(r"^\s*\\begin{([a-zA-Z]*?)}"m, line)) != nothing # environment
            name = m.captures[1]
            xnode = Xnode("environment", Dict("name"=>name))
            rng = findnext(r"^\s*\\begin{([a-zA-Z]*?)}"m, text, pos)
            pos = rng.stop+1
        end
        #rng = findnext(newlineregex, text, pos)
        #rng==nothing && break
        #@show text[rng]
        #pos = rng.stop+1
    end

    return root

end
