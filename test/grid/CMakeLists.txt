# Grid parsing and geometry regression on one representative AutoGrid fixture.
add_test(NAME test_grid_1fkb
         COMMAND test_grid --dpf ${CMAKE_SOURCE_DIR}/data/1fkb/1fkb_epdb.dpf)
