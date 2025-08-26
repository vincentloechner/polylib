#
# $Id: Makefile.am,v 1.3 2007/02/21 21:26:35 skimo Exp $
#

OBJ_DIR = $(top_builddir)

GENERAL =  \
	convex1.in \
	convex2.in \
	gauss1.in \
	herve.in \
	simpl1.in \
	simpl2.in \
	simpl3.in \
	simpl4.in \
	simpl5.in \
	simpl6.in \
	simpl7.in \
	simpl8.in \
	simpl9.in \
	test1.in

EXTRA_DIST = $(GENERAL) \
	$(GENERAL:%.in=%.out) 

tests:
	@failedtest=0; \
	failedtestnames=""; \
	for x in $(GENERAL) ; do \
		file=`basename $$x .in`; \
		printf "Verify file $$file... " ; \
		$(OBJ_DIR)/polytest$(TEST_BITS) < $(srcdir)/$$x > xyz; \
		diff -w xyz $(srcdir)/$$file.out ; \
		result=$$?; \
		if [ "$$result" -eq "0" ]; then \
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
	else \
		echo "Test polytest (general) works correctly"; \
	fi

CLEANFILES=xyz
