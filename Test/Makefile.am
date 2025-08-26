EXTRA_DIST = cases old others
SUBDIRS = Zpolytest general pp ehrhart ehrhart_union ranking

tests:
	@for BITS in $(ALL_BITS); do \
	    printf "\n------------------ Performing tests (bits=$$BITS) ------------------\n"; \
	    for d in $(SUBDIRS); do \
				echo "-- Entering $$d --"; \
				(cd $$d; $(MAKE) TEST_BITS=$$BITS tests) || exit 1;\
	    done \
	done

longtest:
	@for BITS in $(ALL_BITS); do \
	    printf "\n----------- Performing Ehrhart longtests (bits=$$BITS) -----------\n"; \
	    (cd ehrhart; $(MAKE) TEST_BITS=$$BITS longtest) \
	done

dist-hook:
	rm -rf `find $(distdir) -name CVS`
