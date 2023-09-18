include config.mk

.PHONY: all
all: tmp/single-sim_setup single-sim_scripts

.PHONY: single-sim_setup single-sim single-sim_plots single-sim_bifurcation
single-sim_setup single-sim single-sim_plots single-sim_bifurcation: %: tmp/%

tmp/single-sim_setup: src/single-sim_setup.jl
	julia $^
	@touch $@

tmp/single-sim: tmp/single-sim_setup
	julia scripts/single-sim.jl
	@touch $@






ALL_SINGLE_SIM_SCRIPTS = $(wildcard scripts/single-sim_*.jl)
SINGLESIM_SCRIPTS = $(patsubst scripts/%, tmp/%, $(ALL_SINGLE_SIM_SCRIPTS))

$(SINGLESIM_SCRIPTS): tmp/%: scripts/% tmp/single-sim
	julia $<
	@touch $@

.PHONY: single-sim_scripts
single-sim_scripts: $(SINGLESIM_SCRIPTS)




.PHONY: clean
clean:
	rm -rf data/singlesim/*
	rm -rf tmp/*
	rm -rf plots/*
	# rm -rf data/*
