EXTRA_DIST = cases old others
SUBDIRS = Zpolytest ehrhart ehrhart_union ranking general pp

tests:
	@for BITS in $(ALL_BITS); do \
	    echo "\n------------------ Performing tests (bits=$$BITS) ------------------"; \
	    for d in $(SUBDIRS); do \
				(cd $$d; $(MAKE) TEST_BITS=$$BITS tests) || exit 1;\
	    done \
	done

longtest:
	@for BITS in $(ALL_BITS); do \
	    echo "******** Performing tests ($$BITS) ********"; \
	    (cd ehrhart; $(MAKE) TEST_BITS=$$BITS longtest) \
	done

dist-hook:
	rm -rf `find $(distdir) -name CVS`
