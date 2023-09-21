include config.mk

.PHONY: all
all: single-sim-targets ensemble-targets

# Single simulation targets
.PHONY: single-sim_scripts
single-sim_scripts: $(SINGLESIM_SCRIPTS)

SINGLESIM_TARGETS = single-sim_setup single-sim single-sim_plots single-sim_bifurcation
.PHONY: $(SINGLESIM_TARGETS) single-sim-targets
$(SINGLESIM_TARGETS): %: tmp/%
single-sim-targets: $(SINGLESIM_TARGETS)

tmp/single-sim_setup: src/single-sim_setup.jl
	julia $^
	@touch $@

tmp/single-sim: tmp/single-sim_setup
	julia scripts/single-sim.jl
	@touch $@

tmp/single-sim_plots: scripts/single-sim_plots.jl tmp/single-sim
	julia $<
	touch $@

tmp/single-sim_bifurcation: scripts/single-sim_bifurcation.jl tmp/single-sim
	julia $<
	touch $@



# Ensemble targets
ENSEMBLE_TARGETS = ensemble-sim ensemble-sim_single_scenario ensemble-diag-testing_scenarios_plots
.PHONY: $(ENSEMBLE_TARGETS) ensemble-targets
$(ENSEMBLE_TARGETS): %: tmp/%
ensemble-targets: $(ENSEMBLE_TARGETS)

tmp/ensemble-sim: scripts/ensemble-sim.jl
	julia $^
	@touch $@

tmp/ensemble-sim_single-scenario: scripts/ensemble-sim_single-scenario.jl tmp/ensemble-sim
	julia $<
	@touch $@

tmp/ensemble-diag-testing_scenarios_plots: scripts/ensemble-diag-testing_scenarios_plots.jl tmp/ensemble-sim
	julia $<
	@touch $@



.PHONY: clean
clean:
	rm -rf data/singlesim/*
	rm -rf tmp/*
	@echo "cleaning plot output files"
	$(shell fd . 'plots/' -tf | xargs rm -r)
	rm -rf data/*
