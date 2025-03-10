include config.mk

.PHONY: all
all: single-sim-targets ensemble-targets test-targets

# Single simulation targets
SINGLESIM_TARGETS = single-sim_setup single-sim single-sim_plots #single-sim_bifurcation
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
ENSEMBLE_TARGETS = ensemble-sim ensemble-sim_single-scenario ensemble-diag-testing_optimal-thresholds ensemble-diag-testing_constant-thresholds ensemble-diag-testing_scenarios_plots ensemble-diag-testing_optimal-thresholds_single-timeseries
.PHONY: $(ENSEMBLE_TARGETS) ensemble-targets
$(ENSEMBLE_TARGETS): %: tmp/%
ensemble-targets: $(ENSEMBLE_TARGETS)

tmp/ensemble-sim: scripts/ensemble-sim.jl
	julia $^
	@touch $@

tmp/ensemble-diag-testing_scenarios_plots: scripts/ensemble-diag-testing_scenarios_plots.jl tmp/ensemble-sim
	julia $<
	@touch $@

tmp/ensemble-sim_single-scenario: scripts/ensemble-sim_single-scenario.jl tmp/ensemble-sim
	julia $<
	@touch $@

tmp/ensemble-diag-testing_optimal-thresholds: scripts/ensemble-diag-testing_optimal-thresholds.jl tmp/ensemble-sim
	julia $<
	@touch $@

tmp/ensemble-diag-testing_constant-thresholds: scripts/ensemble-diag-testing_constant-thresholds.jl tmp/ensemble-sim
	julia $<
	@touch $@

tmp/ensemble-diag-testing_optimal-thresholds_single-timeseries: scripts/ensemble-diag-testing_optimal-thresholds_single-timeseries.jl tmp/ensemble-sim
	julia $<
	@touch $@


# Manuscript targets
MANUSCRIPT_TARGETS = manuscript
.PHONY: $(MANUSCRIPT_TARGETS) manuscript-targets
$(MANUSCRIPT_TARGETS): %: tmp/%
manuscript-targets: $(MANUSCRIPT_TARGETS)

tmp/manuscript: manuscript/manuscript.typ tmp/ensemble-sim  manuscript/scripts/optimal-thresholds_loading.jl manuscript/scripts/optimal-thresholds_plots.jl manuscript/scripts/optimal-thresholds_tables.jl manuscript/scripts/plotting-setup.jl manuscript/scripts/schematic-plot.jl manuscript/scripts/supplemental_plots.jl manuscript/scripts/supplemental_tables.jl
	@echo "Recomputing schematic, and optimal threshold tables and plots"
	julia manuscript/scripts/optimal-thresholds.jl
	@echo "Compiling manuscript"
	typst compile manuscript/manuscript.typ
	@echo "Compiling supplemental appendix"
	typst compile manuscript/supplemental-appendix.typ
	@touch $@

# Test targets
TEST_TARGETS = runtests
.PHONY: $(TEST_TARGETS) test-targets
$(TEST_TARGETS): %: tmp/%
test-targets: $(TEST_TARGETS)

tmp/runtests: test/runtests.jl single-sim-targets ensemble-targets
	julia $<
	@touch $@

# Cleaning targets
.PHONY: clean-all clean-tmp clean-all-ensemble clean-ensemble-scenarios clean-plots clean-ensemble-sims clean-ensemble-quantiles clean-single-sim clean-ensemble-optimal-thresholds
clean-all: clean-tmp clean-single-sim clean-plots clean-all-ensemble clean-tests

clean-tmp:
	@echo "cleaning all tmp files"
	$(shell rm -rf tmp)

clean-manuscript:
	@echo "cleaning manuscript"
	$(shell fd '.*pdf' './manuscript/' -d 1 -t f -I | xargs rm)
	$(shell fd '.*_files' './manuscript/' -d 2 -t d | xargs rm -r)
	@echo "cleaning manuscript tmp files"
	$(shell fd 'manuscript.*' 'tmp/' | xargs rm)

clean-single-sim:
	@echo "cleaning single-sim output files"
	$(shell fd -g 'single-sim*.jld2' 'data/' -HI | xargs rm -r)
	@echo "cleaning single-sim plot files"
	$(shell fd -g '*.png' 'plots/singlesim/' | xargs rm -r)
	@echo "cleaning single-sim tmp files"
	$(shell fd -g 'single-sim*' 'tmp/' | xargs rm)

clean-plots:
	@echo "cleaning plot output files"
	$(shell fd -g '*.png' 'plots/' -HI | xargs rm -r)

clean-all-ensemble: clean-ensemble-sims clean-ensemble-quantiles clean-ensemble-scenarios clean-ensemble-optimal-thresholds
	@echo "cleaning all ensemble output files"
	$(shell fd . 'out' -td --exclude 'singlesim' -HI | xargs rm -r)
	@echo "cleaning all ensemble tmp files"
	$(shell fd 'ensemble' 'tmp/' | xargs rm -r)

clean-ensemble-sims:
	@echo "cleaning ensemble simulation output files"
	$(shell fd -g 'ensemble-solution*.jld2' 'out/' -HI | xargs rm -r)
	@echo "cleaning ensemble simulation tmp files"
	$(shell fd -g 'ensemble-sim' 'tmp/' | xargs rm)

clean-ensemble-quantiles:
	@echo "cleaning ensemble quantiles"
	$(shell fd -g 'ensemble-quantiles*.jld2' 'out/' -HI | xargs rm -r)
	@echo "cleaning ensemble simulation tmp files"
	$(shell fd -g 'ensemble-sim' 'tmp/' | xargs rm)

clean-ensemble-scenarios:
	@echo "cleaning ensemble scenario files"
	$(shell fd -g 'ensemble-scenario*.jld2' 'out/' -HI | xargs rm -r)
	@echo "cleaning ensemble scenario tmp files"
	$(shell fd -g 'ensemble*scenario' 'tmp/' | xargs rm)
	@echo "cleaning ensemble scenario plot files"
	$(shell fd -g -I '*.png' 'plots/ensemble' | xargs rm -r)

clean-ensemble-optimal-thresholds:
	@echo "cleaning ensemble optimal threshold results"
	$(shell fd -g 'optimal-threshold-results' 'out/' | xargs rm -r)
	@echo "cleaning ensemble optimal threshold plots"
	$(shell fd -g 'optimal-thresholds' 'plots/' | xargs rm -r)
	@echo "cleaning ensemble optimal thresholds tmp files"
	$(shell fd -g 'ensemble-diag-testing_optimal-thresholds' 'tmp/' | xargs rm)

clean-tests:
	@echo "cleaning tests"
	$(shell fd -g 'test' 'tmp/' | xargs rm)
