#
# $Id: Makefile.am,v 1.2 2005/10/20 11:35:25 verdoolaege Exp $
#

OBJ_DIR = $(top_builddir)

PP =  \
	c6.in \
	salvatore.in \
	salvatore_new1.in \
	salvatore_new2.in \
	xavi.in \
	ray1.in \
	ray2.in \
	ray3.in

# this one is pretty long...
#	c7.in \
#

# check those
#	c5.in \
#	c4.in \
#	c9.in
#	c10.in

EXTRA_DIST = $(PP) \
	$(PP:%.in=%.out) 

tests:
	@failedtest=0; \
	failedtestnames=""; \
	for x in $(PP) ; do \
		file=`basename $$x .in`; \
		printf "Verify file $$file... " ; \
		$(OBJ_DIR)/pp$(TEST_BITS) < $(srcdir)/$$x > xyz; \
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
		echo "Test pp works correctly"; \
	fi


CLEANFILES=xyz
