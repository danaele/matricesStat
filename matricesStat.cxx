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

//Calculate nb composante to have 90% of cumulative variance
int numberOfComposantes( std::vector<float> eigenValues)
{
   std::vector<float>::const_iterator eit,eend ;
   float sumEigenValues = 0 ;
   for(eit = eigenValues.begin(), eend = eigenValues.end() ; eit != eend ; eit ++)
   {
       sumEigenValues += *eit;
   }
  std::cout<<sumEigenValues<<std::endl;


  int nbCompo = 0;
  float cumulativeVariance = 0 ;
  for(eit = eigenValues.begin(), eend = eigenValues.end() ; eit != eend ; eit ++)
  {
        nbCompo ++;
        cumulativeVariance  +=  *eit / sumEigenValues ;
        if(sumEigenValues > 0.90)
        {
            return nbCompo;
        }
  }
}

std::vector<float> matrixAsVector ( std::vector< std::vector<float> > matrix )
{
    std::vector<float> matAsVector;
    std::vector< std::vector<float> >::const_iterator it, end;

    for(it=matrix.begin(), end=matrix.end() ; it  !=  end ; it++)
    {
        std::vector<float>::const_iterator vit, vend;
        std::vector<float> row = *it;
        for(vit=row.begin(), vend=row.end() ; vit  !=  vend ; vit++)
        {
            matAsVector.push_back(*vit);
        }
    }
    return matAsVector;
}

int main ( int argc, char *argv[] )
{
  PARSE_ARGS ;
  
  std::ifstream inputFile ;
  std::string labelLine;

  std::vector<std::string> listMatPath;
  std::vector<std::string>::const_iterator lit, lend;

  //list of path for all subjects connectivity matrix
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
  //print_matrix(averageMatrix);
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
  //print_matrix(varianceMatrix);
  write_matrixFile(varianceMatrix,"fdt_network_matrix_variance");


  //---PCA

  //Create each matrix as a vector
  std::list < std::vector<float> > listMatAsVector;
  std::list < std::vector<float> >::const_iterator vit, vend;
  for (it = listMatrix.begin(), end=listMatrix.end() ; it != end ; it++)
  {
      std::vector< std::vector<float> > mat = *it;
      std::vector<float> vectMat;
      vectMat = matrixAsVector(mat);
      listMatAsVector.push_back(vectMat);
  }


  //Matrix of allvectors
  std::vector< std::vector<float> > MatVectors;
  for (vit = listMatAsVector.begin(), vend=listMatAsVector.end() ; vit != vend ; vit++)
  {
      MatVectors.push_back(*vit);
  }
  std::cout<<MatVectors.size()<<std::endl;
  std::cout<<"All matrix as vector"<<std::endl;
  //print_matrix(MatVectors);


  //Create table for PCA with vtk
  float nMat = 0;
  vtkSmartPointer<vtkTable> datasetTable = vtkSmartPointer<vtkTable>::New();
  for (it = listMatrix.begin(), end=listMatrix.end() ; it != end ; it++)
  {
      vtkSmartPointer<vtkDoubleArray> datasetArr = vtkSmartPointer<vtkDoubleArray>::New();

      std::string nameV = "M" + FloatToString(nMat);
      const char* mName = nameV.c_str();
      datasetArr->SetNumberOfComponents(1);
      datasetArr->SetName( mName );
      std::vector< std::vector<float> > mat = *it;
      for(int i= 0 ; i < sizeLine ; i++)
      {
          for(int j= 0 ; j < sizeLine ; j++)
          {

            datasetArr->InsertNextValue(mat.at(i).at(j));
          }
       }
      nMat += 1;
      datasetTable->AddColumn(datasetArr);
      //std::cout<<datasetArr.GetPointer()->GetDataSize()<<std::endl;
  }
   std::cout<<datasetTable.GetPointer()->GetNumberOfColumns() <<std::endl;
   std::cout<<datasetTable.GetPointer()->GetNumberOfRows() <<std::endl;

   int numberValues = datasetTable.GetPointer()->GetNumberOfRows();

  vtkSmartPointer<vtkPCAStatistics> pcaStatistics = vtkSmartPointer<vtkPCAStatistics>::New();
  pcaStatistics->SetInputData( vtkStatisticsAlgorithm::INPUT_DATA, datasetTable );

  for (int i =0 ; i< nMat ; i++)
  {
      std::string nameV = "M" + FloatToString(i);
      const char* mName = nameV.c_str();
      std::cout<<mName<<std::endl;
      pcaStatistics->SetColumnStatus(mName, 1 );
  }
  pcaStatistics->RequestSelectedColumns();
  pcaStatistics->SetDeriveOption(true);
  pcaStatistics->Update();


  //Eigenvalues
   vtkSmartPointer<vtkDoubleArray> eigenvalues = vtkSmartPointer<vtkDoubleArray>::New();
   pcaStatistics->GetEigenvalues(eigenvalues);
   std::vector<float> eigenValues ;
   for(vtkIdType i = 0; i < eigenvalues->GetNumberOfTuples(); i++)
   {
       eigenValues.push_back(eigenvalues->GetValue(i));
       std::cout << "Eigenvalue " << i << " = " << eigenvalues->GetValue(i) << std::endl;
   }



   //Eigenvectors
    vtkSmartPointer<vtkDoubleArray> eigenvectors = vtkSmartPointer<vtkDoubleArray>::New();
    std::vector <std::vector<float> > eigenVectors ;
     pcaStatistics->GetEigenvectors(eigenvectors);
     for(vtkIdType i = 0; i < eigenvectors->GetNumberOfTuples(); i++)
     {
       std::cout << "Eigenvector " << i << " : ";
       double* evec = new double[eigenvectors->GetNumberOfComponents()];
       eigenvectors->GetTuple(i, evec);
       std::vector<float> vector ;
       for(vtkIdType j = 0; j < eigenvectors->GetNumberOfComponents(); j++)
         {
         std::cout << evec[j] << " ";
         vtkSmartPointer<vtkDoubleArray> eigenvectorSingle =
           vtkSmartPointer<vtkDoubleArray>::New();
         pcaStatistics->GetEigenvector(i, eigenvectorSingle);
         vector.push_back(evec[j]);
         }
       delete evec;
       std::cout << std::endl;
       eigenVectors.push_back(vector);
     }

    //Cumulative Variance explained
    int nbCompo = numberOfComposantes(eigenValues);  //90% cumulative variance 
    std::cout<<"Nb compo"<<nbCompo<<std::endl;  //number of eigenVector kept for reconstruction

    Eigen::MatrixXd eigenVector(nbMatrix,nbCompo);
    for(int i = 0; i < nbCompo ; i++)
    {
        for(int j = 0 ; j < nbMatrix ; j++)
        {
            eigenVector(j,i) = eigenVectors.at(i).at(j);
        }
    }

    Eigen::MatrixXd allData(numberValues,nbMatrix);
    for(int i = 0; i < nbMatrix ; i++)
    {
        for(int j = 0 ; j < numberValues ; j++)
        {
            allData(j,i) = MatVectors.at(i).at(j);
        }
    }

    //Reconstruct dataset
    //std::cout<<"Eigenvector size"<<eigenVector.rows()<<"x"<<eigenVector.cols()<<std::endl;
    Eigen::MatrixXd compactdata =eigenVector.transpose() * allData.transpose();
 //  std::cout << "Here is the matrix:\n" << compactdata << std::endl;

 //  std::cout<<compactdata.rows()<<"x"<<compactdata.cols()<<std::endl;

   Eigen::MatrixXd temp =   eigenVector*compactdata;
   Eigen::MatrixXd approx = temp.transpose();

   //std::cout << "Here is the matrix:\n" << approx << std::endl;
   std::cout<<"approx size"<<approx.rows()<<"x"<<approx.cols()<<std::endl;

   //Mean reconstruction imageAllmat
   std::vector < float > vectorMeanMat ;
   for(int i = 0 ; i < approx.rows() ; i++)
   {
       float meanVal = 0;
       for(int j = 0 ; j < approx.cols() ; j++)
       {
           meanVal += approx(i,j);
       }
       meanVal = meanVal / approx.cols();
       vectorMeanMat.push_back(meanVal);
    }

   //Reconstruct matrix
   int nbseed = sqrt(vectorMeanMat.size());
   std::cout<<"nbseed"<<nbseed<<std::endl;

   std::vector < std::vector <float > > matReconstructWithPCA ;

    int val = 0 ;
    for(int i = 0 ; i < nbseed ; i++)
    {
        std::vector <float > line;
        for(int j = 0 ; j < nbseed ; j++)
        {
            line.push_back(vectorMeanMat.at(val));
            val ++;
        }
        matReconstructWithPCA.push_back(line);
     }

     print_matrix(matReconstructWithPCA);
    write_matrixFile(matReconstructWithPCA,"PCAreconstruction");
    
  return 0;

  }

