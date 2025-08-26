#
# $Id: Makefile.am,v 1.2 2005/10/20 11:35:25 verdoolaege Exp $
#

OBJ_DIR = $(top_builddir)

RANKING =  \
	asci03-core.in

EXTRA_DIST = $(RANKING) \
	$(RANKING:%.in=%.out) 

tests:
	@failedtest=0; \
	failedtestnames=""; \
	for x in $(RANKING) ; do \
		file=`basename $$x .in`; \
		printf "Verify file $$file... " ; \
		$(OBJ_DIR)/ehrhart_ranking$(TEST_BITS) < $(srcdir)/$$x > xyz; \
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
		echo "Test ehrhart_ranking works correctly"; \
	fi


CLEANFILES=xyz
