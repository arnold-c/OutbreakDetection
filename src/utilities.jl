export manuscript_dir,
    manuscript_plots,
    supplement_plots,
    supplement_tables


manuscript_dir(args...) = DrWatson.projectdir("manuscript", args...)
manuscript_plots(args...) = manuscript_dir("manuscript_files", "plots", args...)
supplement_plots(args...) = manuscript_dir("supplemental_files", "plots", args...)
supplement_tables(args...) = manuscript_dir("supplemental_files", "tables", args...)
