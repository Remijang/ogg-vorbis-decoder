CXX = g++ CXXFLAGS = -std = c++ 17 - O3 - Wall

                                              all : main

                                                        main : main.cpp
                                                               $(CXX) $(CXXFLAGS) -
                            o main main.cpp

                                clean : rm -
                            f main

                                .PHONY : all clean
