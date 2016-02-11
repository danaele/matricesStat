#include "matricesStat.h"

std::vector< std::vector<float> >read_probtrackx2_matrix( std::string inputMatrixTextFile )
{
    std::string fileName = inputMatrixTextFile;
    std::ifstream inputFile ;
    inputFile.open( fileName.c_str() , std::ios::in ) ;

    std::vector< std::string > mat_string;
    std::vector< std::string >::const_iterator mat_it,mat_end;
    std::string labelLine;

    //Catch each row of the matrix in a vector of string
    if(inputFile.good())
    {
        do
        {
            getline( inputFile , labelLine ) ; //get information line
        }while( labelLine[0] == '#' ) ;

        while(!inputFile.eof())
        {
            if(labelLine.empty() !=true )
            {
                mat_string.push_back(labelLine);
            }
            getline(inputFile,labelLine);
        }
    }
    inputFile.close();

    //Number of row
    int nb_lines=mat_string.size();
    std::cout<<"Number of rows : "<<nb_lines<<std::endl;

    std::vector< std::vector <std::string> > mat;
    std::string::const_iterator s_it, s_end ;

    std::vector< std::vector<float> > mat_num;

    //Find each matrix components
    for(mat_it=mat_string.begin() , mat_end=mat_string.end() ; mat_it!=mat_end ; ++mat_it)
    {
        std::string line = *mat_it;
        std::vector< std::string > cells;
        std::string number;

        for(s_it=line.begin() , s_end=line.end() ; s_it!=s_end ; ++s_it)
        {
            if(*s_it != ' ')
            {
                number += *s_it ;
                if(s_it == s_end-1)  //if last character is a number
                {
                    cells.push_back(number);
                    number.clear();
                }
            }
            else
            {
                if( number != "")
                {
                    cells.push_back(number);
                    number.clear();
                }
                else
                {
                    //Catch next char
                }
            }
        }
        //Verified that is a square matrix
        if(cells.size() != nb_lines)
        {
            std::cout<<"Error : dimension of matrix not matched - Matrix should be square "<<std::endl;
            mat_num.clear();
            return mat_num;
        }
        else
        {
            mat.push_back(cells);

        }
        cells.clear();
    }
    std::cout<<"Matrix dimension verified "<<nb_lines<<"x"<<nb_lines<<std::endl;


    //Convert matrix of string with matrix of double
    for(int i = 0 ; i < nb_lines ; i++)
    {
        std::vector<float>  row_num;
        for(int j = 0 ; j < nb_lines ; j++)
        {
            std::string component = mat[i][j];
            double num_component = atoi( component.c_str() );
            row_num.push_back(num_component);

        }
        mat_num.push_back(row_num);
    }
    return mat_num;
}

void print_matrix(std::vector< std::vector<float> > matrix)
{
    std::cout<<"Matrix : "<<std::endl;
    //Print matrix
    std::vector< std::vector <float> >::const_iterator row_it, row_end;

    for(row_it=matrix.begin() , row_end=matrix.end() ; row_it!=row_end ; ++row_it)
    {
        std::vector< float > row = *row_it;
        std::vector <float>::const_iterator column_it, column_end;
        for(column_it=row.begin() , column_end=row.end() ; column_it!=column_end ; ++column_it)
        {
            float val = *column_it;
            std::cout<<val<<" ";
        }
        std::cout<<""<<std::endl;
    }
}

void write_matrixFile(std::vector< std::vector<float> >  matrix , std::string filename)
{
    std::string filenameOutput=filename;
    std::ofstream outputFile ;
    std::vector< std::vector <float> >::const_iterator row_it, row_end;

    outputFile.open( filenameOutput.c_str() , std::ios::out ) ;
    if( outputFile.good() )
    {
        for(row_it=matrix.begin() , row_end=matrix.end() ; row_it!=row_end ; ++row_it)
        {
            std::vector< float > row = *row_it;
            std::vector <float>::const_iterator column_it, column_end;
            std::string line="";

            for(column_it=row.begin() , column_end=row.end() ; column_it!=column_end ; ++column_it)
            {
                float val = *column_it;
                line += FloatToString(val);
                line += "  ";
            }
            line += "\n";
            outputFile << line ;
        }
    }
    outputFile.close() ;
    std::cout<<"File created : "<<filename<<std::endl;

}

std::string FloatToString ( float number )
{
    std::string result ;
    std::ostringstream convert ;
    convert << number ;
    result = convert.str() ;
    return result ;
}


int main ( int argc, char *argv[] )
{
  PARSE_ARGS ;
  
  std::ifstream inputFile ;
  std::string labelLine;

  std::vector<std::string> listMatPath;
  std::vector<std::string>::const_iterator lit, lend;

  inputFile.open( listMatrixPath.c_str() , std::ios::in );
  if(inputFile.good())
  {
      while(!inputFile.eof())
      {
        getline( inputFile , labelLine ) ;
        //std::cout<<labelLine<<std::endl;
        if(!labelLine.empty())
        {
            listMatPath.push_back(labelLine);
        }
      }

  }
  inputFile.close();

  //ReadMatrixList
  std::list < std::vector< std::vector<float> > > listMatrix;
  std::list < std::vector< std::vector<float> > >::const_iterator it,end;
  int nbMatrix=0;

  for(lit=listMatPath.begin(), lend=listMatPath.end() ; lit!= lend ; lit++)
  {

      std::cout<<*lit<<std::endl;
      std::vector< std::vector<float> > matrix;
      matrix=read_probtrackx2_matrix(*lit);
      listMatrix.push_back(matrix);
      nbMatrix += 1;

  }
    std::cout<<"Number of matrix : "<<nbMatrix<<std::endl;

    int k = 0;
    int sizeLine = 0;
    int sizeRow = 0;
    //Check all matrix have same dimension
    for (it = listMatrix.begin(), end=listMatrix.end() ; it != end ; it++)
    {
        std::vector< std::vector<float> > mat = *it;
        if(k != 0 && sizeLine != mat.size() )
        {
            std::cout<<"all matrix don't have the same dimension - ERROR"<<std::endl;
            return 0;
        }
        sizeLine = mat.size();
        //std::cout<<sizeLine<<std::endl;
        std::vector < std::vector <float> >::iterator rit,rend;
        for (rit = mat.begin(), rend=mat.end() ; rit != rend ; rit++)
        {
            std::vector<float> row = *rit;
            sizeRow = row.size();
            //std::cout<<sizeRow<<std::endl;
            if(sizeRow != sizeLine)
            {
                std::cout<<"all matrix don't have the same dimension - ERROR"<<std::endl;
                return 0;
            }
        }
        k ++ ;
    }
    std::cout<<"All matrix have same dimension : "<<sizeLine<<"x"<<sizeRow<<std::endl;


  //---Average
  std::vector< std::vector<float> > averageMatrix;
  //std::vector< std::vector<float> >::iterator
  int j = 0;
  for (it = listMatrix.begin(), end=listMatrix.end() ; it != end ; it++)
  {
      if(j==0)
      {
          averageMatrix = *it;
      }
      else
      {
          for(int i= 0 ; i < sizeLine ; i++)
          {
              for(int j= 0 ; j < sizeLine ; j++)
              {
                  std::vector< std::vector<float> > mat = *it;
                  averageMatrix.at(i).at(j) =  averageMatrix.at(i).at(j) + mat.at(i).at(j);
              }
          }
      }
      j++;
  }

  for(int i= 0 ; i < sizeLine ; i++)
  {
      for(int j= 0 ; j < sizeLine ; j++)
      {

          averageMatrix.at(i).at(j) = averageMatrix.at(i).at(j) / nbMatrix ;
      }
  }
  print_matrix(averageMatrix);
  write_matrixFile(averageMatrix,"fdt_network_matrix_average");

  //---Variance
  std::vector< std::vector<float> > varianceMatrix;
  for(int i= 0 ; i < sizeLine ; i++)
  {
      std::vector<float>  row;
      for(int j= 0 ; j < sizeLine ; j++)
      {
          row.push_back(0);
      }
      varianceMatrix.push_back(row);
  }

  for (it = listMatrix.begin(), end=listMatrix.end() ; it != end ; it++)
  {

      std::vector< std::vector<float> > mat = *it;
      for(int i= 0 ; i < sizeLine ; i++)
      {
          for(int j= 0 ; j < sizeLine ; j++)
          {
              float val = (mat.at(i).at(j) - averageMatrix.at(i).at(j)) * (mat.at(i).at(j) - averageMatrix.at(i).at(j));
              //std::cout<<val<<std::endl;
              varianceMatrix.at(i).at(j) = varianceMatrix.at(i).at(j) + val;
          }
       }
  }
  for(int i= 0 ; i < sizeLine ; i++)
  {
      for(int j= 0 ; j < sizeLine ; j++)
      {

          varianceMatrix.at(i).at(j) = varianceMatrix.at(i).at(j) / nbMatrix ;
      }
  }
  print_matrix(varianceMatrix);
  write_matrixFile(averageMatrix,"fdt_network_matrix_variance");


  //---PCA

  //Create each matrix as a vector
  std::list < std::vector<float> > listMatAsVector;
  std::list < std::vector<float> >::const_iterator vit, vend;
  for (it = listMatrix.begin(), end=listMatrix.end() ; it != end ; it++)
  {
      std::vector<float> matAsVector;
      std::vector< std::vector<float> > mat = *it;
      for(int i= 0 ; i < sizeLine ; i++)
      {
          for(int j= 0 ; j < sizeLine ; j++)
          {
              matAsVector.push_back(mat.at(i).at(j));
          }
      }
      listMatAsVector.push_back(matAsVector);
  }

  //Matrix of allvectors
  std::vector< std::vector<float> > MatVectors;
  for (vit = listMatAsVector.begin(), vend=listMatAsVector.end() ; vit != vend ; vit++)
  {
      MatVectors.push_back(*vit);
  }
  std::cout<<MatVectors.size()<<std::endl;
  //print_matrix(MatVectors);


  return 0;



  }

