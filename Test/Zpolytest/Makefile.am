#
# $Id: Makefile.am,v 1.2 2005/10/20 11:35:25 verdoolaege Exp $
#

OBJ_DIR = $(top_builddir)

ZPOLYTESTS =  \
	LatDiff1.in \
	LatDiff2.in \
	LatDiff3.in \
	LatDiff3b.in \
	LatDiff5.in \
	LatDiff7.in \
	LatDiff7b.in \
	LatDiff8.in \
	LatDiff8b.in \
	LatHNF1.in \
	LatHNF2.in \
	LatHNF3.in \
	LatInter1.in \
	LatInter2.in \
	LatInter4.in \
	LatInter5.in \
	LatInc1.in \
	LatInc2.in \
	LeftHermite1.in \
	LeftHermite2.in \
	ZEmpty1.in \
	ZEmpty2.in \
	ZInter1.in \
	ZInter2.in  \
	ZInter3.in  \
	ZUnion1.in \
	ZUnion2.in \
	ZDiff1.in \
	ZDiff2.in \
	ZDiff3.in \
	ZDiff4.in \
	ZDiff5.in \
	ZDiff6.in \
	ZDiff7.in \
	ZDiff8.in \
	ZDiff9.in \
	ZImage1.in \
	ZPre1.in \
	ZPre2.in  \
	ZPre3.in \
	ZImPre1.in \
	ZImPre2.in \
	ZImPre3.in

EXTRA_DIST = $(ZPOLYTESTS) \
	$(ZPOLYTESTS:%.in=%.out) 


tests:
	@failedtest=0; \
	failedtestnames=""; \
	for x in $(ZPOLYTESTS) ; do \
		file=`basename $$x .in`; \
		printf "Verify file $$file... " ; \
		$(OBJ_DIR)/Zpolytest$(TEST_BITS) < $(srcdir)/$$x > xyz; \
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
		echo "Test Zpolytest works correctly"; \
	fi

CLEANFILES=xyz
