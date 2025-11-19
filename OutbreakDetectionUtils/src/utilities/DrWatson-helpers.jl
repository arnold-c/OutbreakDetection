using DrWatson: DrWatson

export outdir

outdir(args...) = DrWatson.projectdir("out", args...)
