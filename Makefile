include config.mk

.PHONY: all
all: tmp/single-sim_setup single-sim_scripts ensemble-targets

# Single simulation targets
.PHONY: single-sim_scripts
single-sim_scripts: $(SINGLESIM_SCRIPTS)

SINGLESIM_TARGETS = single-sim_setup single-sim
.PHONY: $(SINGLESIM_TARGETS)
$(SINGLESIM_TARGETS): %: tmp/%


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



# Ensemble targets
ENSEMBLE_TARGETS = ensemble-sim ensemble-diag-testing_scenarios ensemble-sim_single-scenario ensemble-sim_single_scenario_plots ensemble-diag-testing_scenarios_plots
.PHONY: $(ENSEMBLE_TARGETS) ensemble-targets
$(ENSEMBLE_TARGETS): %: tmp/%
ensemble-targets: $(ENSEMBLE_TARGETS)

tmp/ensemble-sim: scripts/ensemble-sim.jl
	julia $^
	@touch $@

tmp/ensemble-diag-testing_scenarios: scripts/ensemble-diag-testing_scenarios.jl tmp/ensemble-sim
	julia $<
	@touch $@

tmp/ensemble-sim_single-scenario: scripts/ensemble-sim_single_scenario.jl tmp/ensemble-sim tmp/ensemble-diag-testing_scenarios
	julia $<
	@touch $@

tmp/ensemble-sim_single_scenario_plots: scripts/ensemble-sim_single_scenario_plots.jl tmp/ensemble-sim_single-scenario
	julia $<
	@touch $@

tmp/ensemble-diag-testing_scenarios_plots: scripts/ensemble-diag-testing_scenarios_plots.jl tmp/ensemble-diag-testing_scenarios
	julia $<
	@touch $@



.PHONY: clean
clean:
	rm -rf data/singlesim/*
	rm -rf tmp/*
	@echo "cleaning plot output files"
	$(fd . 'plots/' -ft | xargs rm)
	# rm -rf data/*
