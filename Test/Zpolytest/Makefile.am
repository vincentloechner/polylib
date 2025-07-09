#
# $Id: Makefile.am,v 1.2 2005/10/20 11:35:25 verdoolaege Exp $
#

OBJ_DIR = $(top_builddir)

ZPOLYTESTS =  \
	LatHNF1.in \
	LatHNF2.in \
	LatHNF3.in \
	LatInter1.in \
	LatInter2.in \
	LatInter3.in \
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
	ZImage1.in \
	ZPre1.in \
	ZPre2.in  \
	ZPre3.in \
	ZImPre1.in \
	SplitZpol1.in \
	SplitZpol2.in \
	SplitZpol3.in

EXTRA_DIST = $(ZPOLYTESTS) \
	$(ZPOLYTESTS:%.in=%.out) 

TEST=$(ZPOLYTESTS:%.in=%.tested) 

tests:
	@failedtest=0; \
	failedtestnames=""; \
	for x in $(ZPOLYTESTS) ; do \
		file=`basename $$x .in`; \
		# echo "$(OBJ_DIR)/Zpolytest$(TEST_BITS) < $(srcdir)/$$x;" \
		echo -n "Verify file $$file... " ; \
		$(OBJ_DIR)/Zpolytest$(TEST_BITS) < $(srcdir)/$$x > xyz; \
		diff -w xyz $(srcdir)/$$file.out ; \
		result=$$?; \
		if [ "$$result" -eq "0" ]; then \
			echo "$$out passed"; \
		else \
			echo "Error: $$file.out is not the same"; \
			failedtest=`expr $$failedtest + 1`; \
			failedtestnames="$$failedtestnames:$$file"; \
	    fi; \
	done ; \
	if [ $$failedtest != 0 ]; then \
		echo "$$failedtest tests failed"; \
		echo "These tests failed: $$failedtestnames:"; \
	else \
		echo "Test Zpolytest works correctly"; \
	fi

CLEANFILES=xyz
