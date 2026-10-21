include(test_macros)
# libgfortran on MinGW restores the numeric locale after each formatted transfer
# from a name the UCRT has freed (GfortranLocale.c). ElmerSolver redirects that
# import at start and says so in its banner; LocaleCheck runs 200 000 formatted
# writes with allocation traffic, after a few percent of which the C runtime's
# numeric locale is the user's without the workaround. A libgfortran that no
# longer imports setlocale fails the banner check: then the workaround can go.
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 line)
RUN_ELMER_EXPECT(EXIT ZERO MATCH_ONCE "libgfortran locale workaround installed"
  MATCH "changed after 0 of 200000 formatted writes" NOMATCH "locale workaround (could not|not installed)")
