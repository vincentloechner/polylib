# Vincent Loechner, 2025

tests:
	@failedtest=0; \
	failedtestnames=""; \
	for x in $(TEST_FILES) ; do \
		file=`basename $$x .in`; \
		printf "Verify file $$file... " ; \
		$(TEST_EXE)$(TEST_BITS) < $(srcdir)/$$x > xyz; \
		if diff -w xyz $(srcdir)/$$file.out ; then \
			printf "passed\n"; \
		else \
			printf "\033[31mError: $$file.out is not the same\033[0m\n"; \
			failedtest=`expr $$failedtest + 1`; \
			failedtestnames="$$failedtestnames:$$file"; \
	    fi; \
	done ; \
	if [ $$failedtest != 0 ]; then \
		echo "\033[31m$$failedtest tests failed\033[0m"; \
		echo "These tests failed: $$failedtestnames:"; \
		exit 1; \
	fi

EXTRA_DIST = $(TEST_FILES) \
	$(TEST_FILES:%.in=%.out) 

CLEANFILES=xyz
