using Plots
using StatsPlots
using DataFrames

f=open("./files/methane_blue.txt")
ch4_blue    = readlines(f)
close(f)

f=open("./files/methane_green.txt")
ch4_green   = readlines(f)
close(f)

f=open("./files/methanol_blue.txt")
ch3oh_blue  = readlines(f)
close(f)

f=open("./files/methanol_green.txt")
ch3oh_green = readlines(f)
close(f)

function get_xy(array)
    df = DataFrame(t=[],x=[])
    for i in array
        s = split(i)
        if isempty(s)
            break
        end
        push!(df, (parse(Float64, s[2]), parse(Float64, s[1])))
    end
    df
end

df1 = get_xy(ch4_blue[23:length(ch4_blue)])
df2 = get_xy(ch4_green[23:length(ch4_green)])
df3 = get_xy(ch3oh_blue[23:length(ch3oh_blue)])
df4 = get_xy(ch3oh_green[23:length(ch3oh_green)])


plot(layout=(2,2))
@df df1 plot!(:x, :t, label="CH4 Front Laser", legend=:bottomright, subplot=1)
@df df2 plot!(:x, :t, label="CH4 Back Laser", color=:green, legend=:bottomright, subplot=2)
@df df3 plot!(:x, :t, label="CH3OH Front Laser",legend=:bottomright, subplot=3)
@df df4 plot!(:x, :t, label="CH3OH Back Laser", color=:green, legend=:bottomright, subplot=4)
title!("Experimental Laser Interferance Pattern", subplot=1)
ylabel!("Laser Intensity (arbitrary units)", subplot=3)
xlabel!("Time (s)", subplot=3)
