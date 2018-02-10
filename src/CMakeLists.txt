ADD_EXECUTABLE(langevin Ronge_Kutta.cpp cubic_spline.cpp dissipative.cpp inertia.cpp langevin.cpp random.cpp shape.cpp store.cpp)
LINK_DIRECTORIES(/usr/local/lib /usr/lib)
TARGET_LINK_LIBRARIES(langevin gsl gslcblas m)