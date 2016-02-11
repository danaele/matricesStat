#include <cstdio> 
#include <matricesStatCLP.h>
#include <itkDirectory.h>
#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <list>
#include <itkImagePCADecompositionCalculator.h>

std::vector< std::vector<float> > read_probtrackx2_matrix( std::string inputMatrixTextFile );
void print_matrix(std::vector< std::vector<float> > matrix);
void write_matrixFile(std::vector< std::vector<float> >  matrix , std::string filename);
std::string FloatToString ( float number );
