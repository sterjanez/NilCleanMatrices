add_library(MatrixLib Matrix.cpp)
target_include_directories(MatrixLib PUBLIC .)

add_library(Matrix8Lib Matrix8.cpp)
target_include_directories(Matrix8Lib PUBLIC .)
target_link_libraries(Matrix8Lib PUBLIC MatrixLib)

add_library(MatrixIntLib MatrixInt.cpp)
target_include_directories(MatrixIntLib PUBLIC .)
target_link_libraries(MatrixIntLib PUBLIC MatrixLib Matrix8Lib)

add_library(PolyMatrixLib PolyMatrix.cpp)
target_include_directories(PolyMatrixLib PUBLIC .)
target_link_libraries(PolyMatrixLib PUBLIC MatrixLib Matrix8Lib MatrixIntLib)

add_library(IntegerSetLib IntegerSet.cpp)
target_include_directories(IntegerSetLib PUBLIC .)
target_link_libraries(IntegerSetLib PUBLIC MatrixLib)

add_library(NilCleanLib NilClean.cpp)
target_include_directories(NilCleanLib PUBLIC .)
target_link_libraries(NilCleanLib PUBLIC MatrixLib Matrix8Lib MatrixIntLib IntegerSetLib)

add_library(UserInterfaceLib UserInterface.cpp)
target_include_directories(UserInterfaceLib PUBLIC .)
target_link_libraries(UserInterfaceLib PUBLIC MatrixLib)

add_library(CommandsLib Commands.cpp)
target_include_directories(CommandsLib PUBLIC .)
target_link_libraries(CommandsLib PUBLIC MatrixLib PolyMatrixLib NilCleanLib UserInterfaceLib IntegerSetLib)

add_executable(Main Main.cpp)
target_include_directories(Main PUBLIC .)
target_link_libraries(Main PUBLIC CommandsLib)

install(TARGETS Main RUNTIME DESTINATION bin)