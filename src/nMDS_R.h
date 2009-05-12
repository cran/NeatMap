#include<iostream>
#include<fstream>
#include <vector>
#include<algorithm>
#include<map>
#include<cmath>
#include<string>

using namespace std;
//String Tokenizer
void Tokenize(const string& str, 
vector<string>& tokens, 
const string& delimiters = "\t") 
{ 
tokens.clear(); 
// Skip delimiters at beginning. 
string::size_type lastPos = str.find_first_not_of(delimiters, 0); 
// Find first "non-delimiter". 
string::size_type pos = str.find_first_of(delimiters, lastPos); 
while (string::npos != pos || string::npos != lastPos) 
{ 
// Found a token, add it to the vector. 
tokens.push_back(str.substr(lastPos, pos - lastPos)); 
// Skip delimiters. Note the "not_of" 
// lastPos = str.find_first_not_of(delimiters, pos); 
lastPos=pos+1; 
if(lastPos==0) 
lastPos=string::npos; 
//lastPos=((pos+1)<string::npos)?pos+1:string::npos; 
// Find next "non-delimiter" 
pos = str.find_first_of(delimiters, lastPos); 
} 
} 
double Pearson(vector <double> data1,vector <double> data2)
{
        double xy=0,x=0,y=0,x2=0,y2=0;
        int N=data1.size();
        
        vector <double>::iterator itr1;
        vector <double>::iterator itr2;
	for(itr1=data1.begin(),itr2=data2.begin();(itr1!=data1.end())&&(itr2!=data2.end());++itr1,++itr2)
	{
		if(isnan(*itr1)||isnan(*itr2))
		{
			N--;
		}
	}

	for(itr1=data1.begin(),itr2=data2.begin();(itr1!=data1.end())&&(itr2!=data2.end());++itr1,++itr2)
    	{
		
		if(!(isnan(*itr1)||isnan(*itr2)))
		{
			xy+=(*itr1)*(*itr2);
			x+=(*itr1);
			y+=(*itr2);
			x2+=(*itr1)*(*itr1);
			y2+=(*itr2)*(*itr2);
     		}
    	}
  
        return -(((xy-x*y/N)/sqrt((x2-x*x/N)*(y2-y*y/N))));
}

double Euclidean(vector <double> data1,vector <double> data2)
{
        double dist=0.0; 
        vector <double>::iterator itr1;
        vector <double>::iterator itr2;
        unsigned int valid=0;

	for(itr1=data1.begin(),itr2=data2.begin();(itr1!=data1.end())&&(itr2!=data2.end());++itr1,++itr2)
    	{
                        if((!isnan(*itr1))&&(!isnan(*itr2)))
                        {
			        dist+=((*itr1)-(*itr2))*((*itr1)-(*itr2));
                                valid++;
                        }
    	}
        return sqrt(dist*valid/data1.size());  
}


unsigned long int remap(unsigned long int n,unsigned long int i,unsigned long int j)
{
  return n*(i)-(i+2)*(i+1)/2+j;
}

class compare
{
public:
double *data;
compare(double* obj) : data(obj) {}

bool operator ()(const unsigned long int p1,const unsigned long int p2)
{
return(data[p1] < data[p2]);
}
};
  


void rank(double *distances,unsigned long int *ranks,unsigned long int size)
{
  vector <unsigned long int> array(size);
  for(unsigned long int i=0;i<size;i++)
  {
    array[ranks[i]]=i;
  }
  sort(array.begin(),array.end(),compare(distances)); 
  for(unsigned long int i=0;i<size;i++)
  {
    ranks[array[i]]=i;
  }
}

class nMDS
{
  
  unsigned long int number_of_points;
  unsigned long int profile_length;
  double *profile_data;
  unsigned long int embed_dim;
  unsigned long int number_of_iterations;
  unsigned int distance_measure;
  
  double *dist_orig;
  double *dist_embed;
  unsigned long int *rank_orig;
  unsigned long int *rank_embed;
  vector <vector <double> > Positions;
  vector <double> Error;

  void Initialize();
  void InitData();
  void CalculateEmbeddedRanks();
  void NormalizePositions();
  void UpdatePositions();

  public:
  
  void Output(double *);
  void ErrorOut(double *);
  
  nMDS(int np,unsigned long int pl,double *data,unsigned long int dim,unsigned int niters,unsigned int dm):number_of_points(np),profile_length(pl),profile_data(data),embed_dim(dim),number_of_iterations(niters),distance_measure(dm)
  {
    Initialize();
    InitData();
    for(unsigned long int iters=0;iters<number_of_iterations;iters++)
    {
        NormalizePositions();
        CalculateEmbeddedRanks();  
        UpdatePositions();
    }
   NormalizePositions();
    
  }

  ~nMDS()
  {
    delete [] dist_orig;
    delete [] dist_embed;
    delete [] rank_orig;
    delete [] rank_embed;
  }
};

void nMDS::Initialize()
{
  dist_orig=new double[number_of_points*(number_of_points-1)/2];
  dist_embed=new double[number_of_points*(number_of_points-1)/2];
  rank_orig=new unsigned long int[number_of_points*(number_of_points-1)/2];
  rank_embed=new unsigned long int[number_of_points*(number_of_points-1)/2];
  vector <double> temp;
  for(unsigned long int i=0;i<number_of_points*(number_of_points-1)/2;i++)
  {
    rank_orig[i]=i;
    rank_embed[i]=i;
    if(i<number_of_points)
    {
        vector <double> temp;
        for(unsigned long int j=0;j<embed_dim;j++)
        {
                temp.push_back(((1.0*random())/RAND_MAX)-0.5);
        }
        Positions.push_back(temp);
    }
  }

}

void nMDS::InitData()
{
        vector <vector <double> > data(number_of_points);
        for(unsigned long int col=0;col<profile_length;col++)
        {
          for(unsigned int row=0;row<number_of_points;row++)
          {
            data[row].push_back(profile_data[row+col*number_of_points]);
          }
        }

        for(unsigned long int i=0;i<number_of_points-1;i++)
        {
          for(unsigned long int j=i+1;j<number_of_points;j++)
          {
                dist_orig[remap(number_of_points,i,j)]=(distance_measure==0)?(Pearson(data[i],data[j])):Euclidean(data[i],data[j]);            
          }
        }
        rank(dist_orig,rank_orig,number_of_points*(number_of_points-1)/2);
}

void nMDS::CalculateEmbeddedRanks()
{
        for(unsigned long int i=0;i<number_of_points-1;i++)
        {
          for(unsigned long int j=i+1;j<number_of_points;j++)
          {
                dist_embed[remap(number_of_points,i,j)]=Euclidean(Positions[i],Positions[j]);            
          }
        }
        rank(dist_embed,rank_embed,number_of_points*(number_of_points-1)/2);
}
void nMDS::NormalizePositions()
{
  vector <double> center_of_mass(embed_dim,0.0);
  double radius=0.0;
  
  for(unsigned long int j=0;j<embed_dim;j++)
  {
        double mean=0.0;
        for(unsigned long int i=0;i<number_of_points;i++)
        {
          mean+=Positions[i][j];
        }
        center_of_mass[j]=mean/number_of_points;
  }
  
  for(unsigned long int i=0;i<number_of_points;i++)
  {
     radius+=Euclidean(center_of_mass,Positions[i]);
  }

  for(unsigned long int i=0;i<number_of_points;i++)
  {
        for(unsigned long int j=0;j<embed_dim;j++)
        {
                Positions[i][j]=(Positions[i][j]-center_of_mass[j])/radius;
        }
  }
        
}
void nMDS::UpdatePositions()
{
        vector <vector <double> > temp=Positions;
        double error=0.0;
        for(unsigned long int i=0;i<number_of_points-1;i++)
        {
          for(unsigned long int j=i+1;j<number_of_points;j++)
          {
                double c=(1.0*rank_orig[remap(number_of_points,i,j)]-1.0*rank_embed[remap(number_of_points,i,j)]);
                error+=abs(c);
                for(unsigned long int d=0;d<embed_dim;d++)
                {
                  temp[i][d]+=(c)*(Positions[i][d]-Positions[j][d])/(number_of_points*number_of_points*number_of_points);
                  temp[j][d]+=(c)*(Positions[j][d]-Positions[i][d])/(number_of_points*number_of_points*number_of_points);
                }
          }
        }
        Positions=temp; 
        Error.push_back(error);
}

void nMDS::Output(double *positions_out)
{
  for(unsigned long int d=0;d<embed_dim;d++)
  {
        for(unsigned long int i=0;i<number_of_points;i++)
        {
          positions_out[i+d*number_of_points]=Positions[i][d];
        }
  }
}

void nMDS::ErrorOut(double *error_out)
{
        for(unsigned int i=0;i<number_of_iterations;i++)
        {
          error_out[i]=Error[i];
        }
}

