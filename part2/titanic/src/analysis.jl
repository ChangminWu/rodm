#using Plots
using CSV
using DataFrames
using Statistics
using StatsPlots, Interact
using Blink

dataSet = "titanic"
dataFolder = "/home/changmin/Documents/Projects/Dascim/rodm/part2/titanic/data/"
file = dataFolder * dataSet * ".csv"

data = CSV.read(file,  header=true)
describe(data)

w = Window()
body!(w, dataviewer(data))
