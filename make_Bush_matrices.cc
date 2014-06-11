#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdio.h>      /* printf */
#include <math.h> 
using namespace std;




//Global variable
//integers
const int s=8;//size of our matrix in the form of 2^k
const int s2 = s*s;
const int k=3;//k

//matrix

int product[s][s];//multiplication table will save in it if we call the function multiplication_table
int H[s][s];//this is a hadamard matrix of size s in the form we use in our report. 
int R[s][s][s];//Made by reapition of rows of generalized Hadamard matrix
int IR[s][s][s]={0};// this is the matrix of form IxIxI, IxIxR, .... we should call the function generat_IR to fill this function bace in index i
int M[s][s][s]={{0}};// the componentwise_multiplication of H and R
//following four variable are for making the Hadamardf matrix
int H1[s][s] = {{0}};
int H2[2*s][2*s] = {{0}};
const int h1[2][2] = {{1,1},{1,-1}};
const int h2[4][4] = {{1,1,1,1},{1,1,-1,-1},{1,-1,1,-1},{1,-1,-1,1}};;
//the following variable used to change the 1234... matrix to 0 1 -1 matrix and change matrix IR to replac. 
int replac[s][s2][s2];


//declaration of function
bool check_weighing(int );//check if the matrix is weighning
bool check_BUSH(int** );//check if the matrix is bush type or not;
void multiplication_table();//save the multiplication table in matrix product using multiplication function
int multiplication(int ,int);//return the multiplication of two arguments using binary code
int*  binary(int );//change the argument to binary
int integer(int* );//recieve a binary and change it to an integer

//the six following functions do the normal matrix multiplication 
int** matrix_multiplication(int[s][s] ,int [s][s]);
int** matrix_multiplication(int[s][s][s] ,int [s][s][s]);
int** matrix_multiplication(int**,int [s][s]);
int** matrix_multiplication(int[s][s] ,int**);
int** matrix_multiplication(int a[s][s2][s2],int** b, int y);
int** matrix_multiplication(int** ,int**);
////return the transposition of the argument
int** Trans(int a[s][s]);
int** Trans(int a[s][s2][s2], int);
//the four following function do kronecker product
int** Kronecker(int** a, int ** b, int i=s, int j=s );
int** Kronecker(int a[s][s], int ** b ,  int j =s);
int** Kronecker(int** a,  int b[s][s] ,  int i =s);
int** Kronecker(int a[s][s], int b[s][s]);

//the infinity many number  following functions print a matrix
void print(const int[s][s]);
void print(int[s][s][s]);
void print(int[s][s2][s2]);
void print(int b[s][s2][s2], int y);
void print1(int[2*s][2*s]);
void print(int**, int i =s);
//
void generat_IR (int,int ,int ,int ,int );//save the IxIxR in IR matrix, the third int shoul be 1.
void generat_HH();//save Hadamard matrix in H1 and H2


void fill_row();//this function fill the variable row
void generat_GH();//this function save a generalized hadamard matrix in variable p[][] which is anarray object of class polynomail
void make_matrix_M();//fill w by componentwise product of H and  R
void make_matrix_R();//this function fill R[s][s][s]
void print_matrix_P();//this matrix print the array P using specific method
void print_row();//print the element of row using specific method
//
void replaced();//this function fill replac varaible
int* generate_all_rows(int index);//This function generates a row of 0,1,-1 depending on the index(varying between 0 and 3^s2) input
bool check_unbiased(int index, int a[s][s2][s2], int a_index, int row_index);//This function checks if the row obtained from the index, and the (row_index)th row of the (a_index)th weighing matrix
// are unbiased.

bool check_symmetric(int** p, int n); //checks if the matrix of order n is symmetric.

//class of polynomail
class polynomial
{
	public:
	int coeff[k];
	int index;
	void compute_index()
	{
		int c=1;
		for(int i=0;i<k;i++)
		{
			c = c + coeff[i]*pow(2,i);
		}
	        index = c;
	}
};


//the two following variables are used to make a generqalized hadamard matrix
polynomial p[s][s];
polynomial row[s];

int main()
{ 
   
   //initilization of h2;
   for(int i= 0; i<4; i++)
   {
      //  h2[0][i]=1;
      //   h2[i][0] = 1;
   }
   // h1[0][0]=1;
   //  h1[0][1]=1;
   //  h1[1][0]=1;
   //  h1[1][1]=-1;

   //  h2[1][1] =1;
   //  h2[1][2] =-1;
   //  h2[1][3] =-1;
//   h2[2][1] =-1;
//   h2[2][2] =1;
//   h2[2][3] =-1;
//   h2[3][1] =-1;
//   h2[3][2] =-1;
//   h2[3][3]=1;
   //
   //for(int i=0;i<s;i++)
   //  {
   //   for(int j=0;j<s;j++)
//	 cout<<[i][j]<<" ";
   //   cout<<endl;
//   }
   //multiplication_table();
   /*  for(int i = 0; i < s; i++)
       {
       for(int j = 0;  j < s; j++)
       {
       cout<<product[i][j]<<" ";
       }
      cout<<endl;
      }*/
   //  int d[s][s]={{0,1,1,1,1,1,1,1},{1,0,1,1,1,1,1,1},{1,1,1,1,1,0,1,1},{1,1,0,1,1,1,1,1},{1,1,1,0,1,1,1,1},{1,1,1,1,0,1,1,1},{1,1,1,1,1,1,0,1},{1,1,1,1,1,1,1,0}} ;
   //  int h[s][s]={{2,0,2,2,2,2,2,9},{2,8,2,2,2,2,4,9},{2,2,6,2,2,2,5,9},{2,2,2,2,2,2,2,9},{2,2,2,2,0,2,4,9},{2,2,2,2,2,2,2,9},{2,9,2,2,2,2,2,9},{1,1,1,1,1,1,1,9} } ;
   // print( Kronecker(d,h), s2);
   
   
/* print( matrix_multiplication(h,d));
   cout<<endl;
   print( matrix_multiplication(Trans(h),d));
   cout<<endl;
   print(matrix_multiplication(h, Trans(d)));
   cout<<endl;
   print(matrix_multiplication(Trans(h),Trans(d)));
*/
   ///// DO NOY change any thing in the foloowing part
   for(int i=0; i < s; i++)
   {      
     IR[i][0][i]=1;
     generat_IR(0,i,1,i,i);
   }
   generat_HH();
   fill_row();
   generat_GH();
   make_matrix_R();
   make_matrix_M();
   replaced();
//////
   /* all 8 matrices of order 64 t0 64 in replac are weghing matrix W(64,8);
   for(int u = 0; u < s; u++)
   {
      cout<<"check "<<u+1<<" :"<<endl;
      check_weighing(u);
      }*/
   ///////
 
   // check_BUSH( matrix_multiplication(replac,Trans(replac,1),0));

  //  print(replac,1);
      //  cout<<endl;
   //  print(Trans(replac,1),64);
//   cout<<endl;
   //  cout<<"Weighing matrix of order 64"<<endl;
   //  print(replac,1);
   //  cout<<"Bush type Hadamard matrix of order 64"<<endl;
   
//   print(matrix_multiplication(replac,Trans(replac,1),0),64);
   // cout<<"Printing h1"<<endl;
   print(replac);
   //print_row();
   // print_matrix_P();
   // print(h2);
  
   

   for(int i=0;i<s;i++)
   {
	   for(int j=i+1;j<s;j++)
	   {
		   if(check_symmetric(matrix_multiplication(replac,Trans(replac,j),i),s2) && check_BUSH(matrix_multiplication(replac,Trans(replac,j),i)))
		   {
			   cout<<"Symmetric Bush type matrix"<<endl;
		   }
		   else
			   cout<<"Un-symmetric Bush type matrix"<<endl;
		   print(matrix_multiplication(replac,Trans(replac,j),i),s2);
	   }
   }

			  
   
   /*
   //The following code checks if a random row of 1,0,-1 generated is unbiased with any of the rows of any of the s weighing matrices.
   long long int index_max,ind;
   int num=0;
   ind = s2;
   cout<<"The current index is"<<endl;
   index_max = pow(3,ind);
   cout<<"Index_max"<<index_max<<endl;
   for(int l=0;l<s;l++)
   {// l is the index for the weighing matrix.
	   for(int m =0;m<s2;m++)
	   {//m is the index for the rows of a weighing matrix.
		   num=0;
		   for(int i=0;i<s2;i++)
		   {
			   cout<<replac[l][m][i]<<" ";
		   }
		   cout<<endl;
		   for(int p=0;p<index_max;p++)
			   {// p is the index of the row being generated, whose unbiasedness with the rows of the existing weighing matrices must be checked.
					   if(check_unbiased(p,replac,l,m)==true)
						{num=num+1;}
			   }
                  //num calculates the number of such unbiased rows found.
		  cout<<"num ="<<num<<endl;
	   }
   }
   */
   return 0;

}

bool check_symmetric(int** p, int n)
{
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			if(p[i][j]!=p[j][i])
				return false;
		}
	}
	return true;
}

int* generate_all_rows(int ind)
{
        int* row = new int[s2];
        int n=ind;
        int digit=0;
        int sum=0;
        int i=0;
        while(i<s2)
        {
                digit = n%3;
                if(digit==2)
                        digit=-1;
                n=n/3;
                row[i] = digit;
                i++;
        }
        return(row);
}

bool check_if_rows_match(int p[s2])
{
	
	for(int i=0;i<s;i++)
	{
		for(int j=0;j<s2;j++)
		{
			for(int k=0;k<s2;k++)
			{
				if(p[i]!=replac[i][j][k])
					return false;
			}
		}
	}
	return true;
}


float magnitude( int p[s2])
{
	float sum=0;
	for(int i=0;i<s2;i++)
	{
		sum=sum + p[i]*p[i];
	}
	return(sqrt(sum));
}
void print(int p[s2])
{
	for(int i=0;i<s2;i++)
	{
		cout<<p[i]<<" ";
	}
	cout<<endl;
}
bool check_unbiased(int ind, int a[s][s2][s2], int a_index, int row_index)
{
        int* row_gen = generate_all_rows(ind);
        float sum=0;
        for(int i=0;i<s2;i++)
        {       
                sum= sum + row_gen[i]*a[a_index][row_index][i];
        }
        float c= magnitude(row_gen);
	cout<<"Magnitude of row_gen is ="<<c<<endl;
	if(c!=0)
	{
		cout<<"The inner product of row from weighing matrix and the row_gen is "<<sum<<endl;
		sum=sum/(sqrt(s)*c);
		cout<<"And the magnitude of the inner product of the unit vectors is"<<sum<<endl;
	}

	if(sum==0 || sum==1 ||sum==-1)
	{
            if(check_if_rows_match(row_gen)==false)
	    {
		    print(row_gen);
		    cout<<"sum ="<<sum<<endl;
		    return true;
	    }
	    else
		    return false;
	}
        else
	{

                return false;
	}
}


bool check_weighing(int d)
{
   for(int i = 0; i < s2; i++)
   {
      for(int j = 0;  j < s2; j++)
      {
	 if(replac[d][i][j] != 1 && replac[d][i][j] != -1 && replac[d][i][j] != 0)
	 {
	  
	    return false;
	 }

	    
      }
   }
   int sum;
   for(int i= 0; i < s2; i++)
   {
      sum = 0;
      for(int j =0; j < s2; j++)
      {
	 sum =sum+ (replac[d][i][j]* replac[d][i][j]);
      }
 
      if(sum != 8)
	 return false;
   }
   for(int i = 0; i < s2; i++)
   {
      for(int j = i+1; j<s2; j++)
      {
	 sum = 0;
	 for(int u = 0; u < s2; u++)
	 {	 
	    sum =sum + (replac[d][i][u]* replac[d][j][u]);
	 }
	 if(sum != 0)
	    return false;
      }
   }
   return true;
}

bool check_BUSH(int** A)
{
   int sum = 0;
   for(int i = 0; i < s2; i++)
   {
      for(int j = 0;  j < s2; j++)
      {
	 if(A[i][j] != 1 && A[i][j] != -1)
	 {
	  
	    return false;
	 }
    
      }
   }
     
   for(int i = 0; i < s2; i++)
   {
      for(int j = i+1; j<s2; j++)
      {
	 sum = 0;
	 for(int u = 0; u < s2; u++)
	 {	 
	    sum =sum + (A[i][u]* A[j][u]);
	 }
	 if(sum != 0)
	    return false;
      }
   }
   //two outer for loop are block while the two inner for are elements of each block
   for(int i1 = 0; i1 < s; i1++)
   {
      for(int i2 = 0; i2 < s; i2++)
      {
	 if(i1 ==i2)
	 {
	    for(int i3 = 0; i3 < s; i3++)
	    {
	       for(int i4 =0; i4 <s; i4++)
	       {
		  if( A[i3+(i1*s)][i4+(i2*s)] != 1)
		     return false;
		  
	       }
       	    }
	 }
	 else
	 {
	    
	    for(int i3 = 0; i3 < s; i3++)
	    {
	       sum =0;
	       for(int i4 =0; i4 <s; i4++)
	       {
		  sum = sum  +  A[i3+(i1*s)][i4+(i2*s)];
	       }
	       if(sum != 0)
		  return false;
       	    }
	 }
      }
   }  
   return true;
}

void multiplication_table()
{
   for(int i= 0; i < s; i++)
   {
      for(int j = 0; j < s; j++)
      {
	 product[i][j] = multiplication(i,j);
      }
   }
 }
int multiplication(int y, int z)
{
   int result = 0;
   int count = 0; 
   int* bin1;
   int* bin2;
   int* bin3=new int[k];
   bin1 = binary(y);
   bin2 = binary(z);
   
     for(int i =0; i < k; i++)
     {
	bin3[i] = (bin1[i]+bin2[i])%2;
     }
     
     return(integer(bin3));
}
int*  binary(int a)
{
  
   int* p = new int[k];
   int count=0;
   while(a!=0)
   {
      p[count] = a%2;
      a = a/2;
      count ++;
   }
   return(p);
}

int integer(int* a)
{
  
   int result=0;
   for(int i= 0; i< k; i++)
     {
	result = result + (a[i]*pow(2.0,i));

     }
    return result;
}


int** matrix_multiplication(int a[s][s],int b[s][s])
{
   
   int** c =new int*[s2];
   for(int i=0; i < s2; i++)
   {
      c[i] = new int[s2];
   }
   int sum;
   for(int i = 0; i < s2; i++)
   {
      for(int j = 0; j < s2; j++)
      {
	 sum =0;
	 for(int o =0; o < s2; o++)
	 {
	    sum = sum + (a[i][o]*b[o][j]);
	 }

	  c[i][j]=sum;
      }
   }
   return (c);
}
int** matrix_multiplication(int a[s][s2][s2] ,int b[s][s2][s2], int y)
{
    int** c =new int*[s2];
   for(int i=0; i < s2; i++)
   {
      c[i] = new int[s2];
   }
   int sum;
   for(int i = 0; i < s2; i++)
   {
      for(int j = 0; j < s2; j++)
      {
	 sum =0;
	 for(int o =0; o < s2; o++)
	 {
	    sum = sum + (a[y][i][o]*b[y][o][j]);
	 }

	  c[i][j]=sum;
      }
   }
   return (c);
}

int** matrix_multiplication(int a[s][s2][s2],int** b, int y)
{
   
   int** c =new int*[s2];
   for(int i=0; i < s2; i++)
   {
      c[i] = new int[s2];
   }
   int sum;
   for(int i = 0; i < s2; i++)
   {
      for(int j = 0; j < s2; j++)
      {
	 sum =0;
	 for(int o =0; o < s2; o++)
	 {
	    sum = sum + (a[y][i][o]*b[o][j]);
	 }

	  c[i][j]=sum;
      }
   }
   return (c);
}
int** matrix_multiplication(int** a,int b[s][s])
{
   
   int** c =new int*[s];
   for(int i=0; i < s; i++)
   {
      c[i] = new int[s];
   }
   int sum;
   for(int i = 0; i < s; i++)
   {
      for(int j = 0; j < s; j++)
      {
	 sum =0;
	 for(int o =0; o < s; o++)
	 {
	    sum = sum + (a[i][o]*b[o][j]);
	 }

	  c[i][j]=sum;
      }
   }
    return (c);
}
int** matrix_multiplication(int** a,int** b)
{
   
   int** c =new int*[s];
   for(int i=0; i < s; i++)
   {
      c[i] = new int[s];
   }
   int sum;
   for(int i = 0; i < s; i++)
   {
      for(int j = 0; j < s; j++)
      {
	 sum =0;
	 for(int o =0; o < s; o++)
	 {
	    sum = sum + (a[i][o]*b[o][j]);
	 }

	  c[i][j]=sum;
      }
   }
  
   return (c);
}

int** Trans(int a[s][s])
{
   int** c =new int*[s];
   for(int i=0; i < s; i++)
   {
      c[i] = new int[s];
   }
   for(int i = 0; i < s; i++)
   {
      for(int j = 0; j < s; j++)
      {
	 c[i][j]= a[j][i];
      }
   }
 
   return c;
}


int** Trans(int a[s][s2][s2], int y)
{
    int** c =new int*[s2];
   for(int i=0; i < s2; i++)
   {
      c[i] = new int[s2];
   }
   for(int i = 0; i < s2; i++)
   {
      for(int j = 0; j < s2; j++)
      {
	 c[i][j]= a[y][j][i];
      }
   }
 
   return c;
}

int** Kronecker(int** a, int ** b, int m, int n)
{
   int q =0;
   int r =0;
 
   int** c =new int*[m*n];
   for(int i=0; i < m*n; i++)
   {
      c[i] = new int[s2];
   }

   for(int i= 0; i < m*n; i++)
   {
      for(int j= 0; j <m*n; j++)
      {
	 c[i][j] = a[i/m][j/m]*b[i%n][j%n];
      }
   }
   
   return c;
}

int** Kronecker(int a[s][s] ,int ** b, int n)
{
 
 
   int** c =new int*[s*n];
   for(int i=0; i < s*n; i++)
   {
      c[i] = new int[s2];
   }

   for(int i= 0; i < s*n; i++)
   {
      for(int j= 0; j <s*n; j++)
      {
	 c[i][j] = a[i/s][j/s]*b[i%n][j%n];
      }
   }
   
   return c;
}

int** Kronecker(int** a, int b[s][s],int m)
{
   int q =0;
   int r =0;
 
   int** c =new int*[s*m];
   for(int i=0; i < s*m; i++)
   {
      c[i] = new int[s2];
   }

   for(int i= 0; i < s*m; i++)
   {
      for(int j= 0; j <s*m; j++)
      {
	 c[i][j] = a[i/m][j/m]*b[i%s][j%s];
      }
   }
   
   return c;
   
}

int** Kronecker(int a[s][s],int b[s][s])
{
   int q =0;
   int r =0;
 
   int** c =new int*[s2];
   for(int i=0; i < s2; i++)
   {
      c[i] = new int[s2];
   }

   for(int i= 0; i < s2; i++)
   {
      for(int j= 0; j <s2; j++)
      {
	 c[i][j] = a[i/s][j/s]*b[i%s][j%s];
      }
   }
   
   return c;
}

void print(const int a[s][s])
{
     for(int i=0;i<s;i++)
   {
      for(int j=0;j<s;j++)
	 cout<<a[i][j]<<" ";
      cout<<endl;
   }
}
void print(int b[s][s][s])
{
   for(int l=0; l < s; l++)
   {
   for(int i=0;i<s;i++)
   {
      for(int j=0;j<s;j++)
	 cout<<b[l][i][j]<<" ";
      cout<<endl;
   }
   cout<<endl;
   cout<<"______________"<<endl;
   }
}

void print(int b[s][s2][s2])
{
   for(int l=0; l <s ; l++)
   {
   for(int i=0;i<s2;i++)
   {
      for(int j=0;j<s2;j++)
	 cout<<b[l][i][j]<<" ";
      cout<<endl;
   }
 
   cout<<"_________________________________________________________________________"<<endl;
     cout<<endl;
   }
}
void print(int b[s][s2][s2], int y)
{
 
   for(int i=0;i<s2;i++)
   {
      for(int j=0;j<s2;j++)
	 cout<<b[y][i][j]<<" ";
      cout<<endl;
   }
 

}

void print1(int a[2*s][2*s])
{
     for(int i=0;i<2*s;i++)
   {
      for(int j=0;j<2*s;j++)
	 cout<<a[i][j]<<" ";
      cout<<endl;
   }
}
void print(int** a , int k)
{
    for(int i=0;i<k;i++)
   {
      for(int j=0;j<k;j++)
	 cout<<a[i][j]<<" ";
      cout<<endl;
   }
}

void generat_IR(int indexi,int indexj,int size,int i,int index_IR)
{
        if(size==s)
                return;
        int q = i%2;
        i=i/2;
        if(q==0)
        {
                for(int m=0;m<size;m++)
                        for(int l=0;l<size;l++)
                        {
			   IR[index_IR][indexi+size+m][indexj+size+l]=IR[index_IR][indexi+m][indexj+l];
			 
                        }
		generat_IR(indexi,indexj,2*size,i,index_IR);
        }
        else if(q==1)
        {
                for(int m=0;m<size;m++)
                        for(int l=0;l<size;l++)
                        {
			   IR[index_IR][indexi+size+m][indexj-size+l]=IR[index_IR][indexi+m][indexj+l];
		
                        }
		generat_IR(indexi,indexj-size,2*size,i,index_IR);
        }
}

void generat_HH()
{
     
   
    for(int i= 0; i < s; i++)
   {
      for(int j= 0; j <s; j++)
      {
	 if(i>=4&&j>=4)
	    H1[i][j] =-h2[i%4][j%4];
	 else
	      H1[i][j] = h2[i%4][j%4];
      }
   }
    
     
    for(int i= 0; i < 16; i++)
   {
      for(int j= 0; j <16; j++)
      {
	 if(i>=8&&j>=8)
	    H2[i][j] =-H1[i%4][j%4];
	 else
	      H2[i][j] = H1[i%4][j%4];
      }
   }
  
    
}

void fill_row()
{
	int n;
	bool binary;
	int i;
	for(int ind=0;ind<s;ind++)
	{
		i=0;
		n=ind;
		while(i<=k)
		{
		   binary = n%2;
		   n=n/2;
		   row[ind].coeff[i]=binary;
		   i++;
		}
	}
}




void generat_GH()
{
   int c;
   if(s!=2)
   {
   for(int ind_i=0;ind_i<s;ind_i++)
   {
      for(int ind_j=0;ind_j<s;ind_j++)
      {
	 for(int i=0;i<k;i++)
	 {
	    for(int j=0;j<k;j++)
	    {
	       c = row[ind_i].coeff[i]*row[ind_j].coeff[j];
	       if((i+j)<k)
		  p[ind_i][ind_j].coeff[i+j]=(p[ind_i][ind_j].coeff[i+j]+c)%2;
	       else
	       {
		  p[ind_i][ind_j].coeff[(i+j)-k+1] = (p[ind_i][ind_j].coeff[(i+j)-k+1] + c)%2;
		  p[ind_i][ind_j].coeff[(i+j)-k] = (p[ind_i][ind_j].coeff[(i+j)-k] +c)%2;
	       }
	    }
	 }
	 p[ind_i][ind_j].compute_index();
      }
   }
   }
   else
   {
      p[0][0].index=1;
      p[0][1].index=1;
      p[1][0].index=1;
      p[1][1].index=2;
   }
   
}

void make_matrix_M()
{
   if(s==2)
  {
    for(int y =0 ; y <s ; y++)
      {
         for(int j=0;j<s;j++)
         {
            for(int l=0;l<s;l++)
            {
             
	       M[y][j][l] = h1[j][l]*R[y][j][l];
            }
	    cout<<endl;
         }
      }
  }
if(s==4)
{
for(int y =0 ; y <s ; y++)
      {
         for(int j=0;j<s;j++)
         {
            for(int l=0;l<s;l++)
            {
               M[y][j][l] = h2[j][l]*R[y][j][l];
            }
         }
      }
}
   if(s==8)
   {
      for(int y =0 ; y <s ; y++)
      {
	 for(int j=0;j<s;j++)
	 {
	    for(int k=0;k<s;k++)
	    {
	       M[y][j][k] = H1[j][k]*R[y][j][k];
	    }
	 }
      }
   }
   if(s==16)
   {
      for(int y =0 ; y <s ; y++)
      {
	 for(int j=0;j<s;j++)
	 {
	    for(int k=0;k<s;k++)
	    {
	       M[y][j][k] = H2[j][k]*R[y][j][k];
	    }
	 }
      }
   }
}

void make_matrix_R()
{
   for(int y =0 ; y <s ; y++)
   {
        for(int j=0;j<s;j++)
        {
                for(int k=0;k<s;k++)
                {
                        R[y][j][k] = p[y][k].index;
                }
        }
   }
}
void print_matrix_P()
{
        cout<<endl;
        cout<<"_________________________________"<<endl;
        for(int ind_i=0;ind_i<s;ind_i++)
        {
                for(int ind_j=0;ind_j<s;ind_j++)
                        {
                                        
                        cout<<p[ind_i][ind_j].index<<" ";
 
                        }
                cout<<endl;
        }

}
void print_row()
{
	for(int i=0;i<s;i++)
	{
		for(int j=0;j<k;j++)
		{
		 	cout<<row[i].coeff[j]<<"T^"<<j<<" + ";
		}
		cout<<"|";
	}
}

void replaced()
{
   for(int i = 0; i < s ; i++)
   {
     
      for(int j= 0; j < s; j++)
      {
	 
	 for(int m =0; m < s; m++)
	 {
	  
	    for(int u= 0; u < s; u++)
	    {
	     
	       for(int z =0; z <s; z++)
	       {
		  if(s==2)
		  {
		     replac[i][u+(s*j)][z+(s*m)] = IR[R[i][j][m]-1][u][z]*h1[j][m];
		  }
		  else if(s==4)
		  {
		     replac[i][u+(s*j)][z+(s*m)] = IR[R[i][j][m]-1][u][z]*h2[j][m];
		  }
		 else if(s==8)
		  {
		     replac[i][u+(s*j)][z+(s*m)] = IR[R[i][j][m]-1][u][z]*H1[j][m];
		  }
		 else if(s==16)
		  {
		     replac[i][u+(s*j)][z+(s*m)] = IR[R[i][j][m]-1][u][z]*H2[j][m];
		     
		  }
	       }
	       
	    }
	 }
      }
   }
}
