# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign subdir-objects

# Don't loose the autoconf include path
ACLOCAL_AMFLAGS = -I config

# Subdirectories to build
SUBDIRS = src 

# Utility target: build but don't run tests
build-check:
	$(MAKE) 'TESTS_ENVIRONMENT=: ' check

doc_DATA = README.md doc/GreatCircleArcIntersections.mw doc/GreatCircleLatitudeIntersections.mw

nobase_dist_bin_SCRIPTS = \
	test/areaint.m \
	test/run.sh \
	test/run_fvtogll_cs.sh \
	test/run_fvtogll_ico_diffnorms.sh \
	test/run_glltofv_rll_diffnorms.sh \
	test/testdata_diffnorms_mono.sh \
	test/generate_testdata.sh \
	test/run_fvtofv_ico.sh \
	test/run_fvtogll_cs_diffnorms.sh \
	test/run_fvtogll_rll.sh \
	test/run_glltofv_rll_m.sh \
	test/optimtest.m \
	test/run_fvtofv_ico_diffnorms.sh \
	test/run_fvtogll_cs_m.sh \
	test/run_fvtogll_rll_diffnorms.sh \
	test/run_glltofv_rll_m_diffnorms.sh \
	test/plot_meshes.m \
	test/run_fvtofv_rll.sh \
	test/run_fvtogll_cs_m_diffnorms.sh \
	test/run_glltofv_ico.sh \
	test/run_glltogll_cs.sh \
	test/plot_nodes.m \
	test/run_fvtofv_rll_diffnorms.sh \
	test/run_fvtogll_ico.sh \
	test/run_glltofv_rll.sh \
	test/run_glltogll_cs_diffnorms.sh

EXTRA_DIST = $(doc_DATA)

CLEANFILES = configs.sed

