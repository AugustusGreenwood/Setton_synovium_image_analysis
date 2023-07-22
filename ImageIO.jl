function get_spacing(path)::Float64
    spacing = ""
    open(path, "r") do file
        for _ in 1:50
            line = replace(readline(file), "\0"=>"")
            name, value = split_around_equals(line)
            if name == "spacing"
                spacing = parse(Float64, value)
            elseif name == "unit"
                if value != "micron"
                    println("WARNING: units are not in microns")
                end 
            end
        end
    end
    return spacing
end


function split_around_equals(str::String)::Tuple{String, String}
    index = collect(findfirst("=", str))[1]
    name = str[1:index-1]
    value = str[index+1:end]
    return name, value
end