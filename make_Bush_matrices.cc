#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdio.h>      /* printf */
#include <math.h> 
using namespace std;


//Things to be done.
//D1. Initialize s=16, k=4, s2=256
//D2. Create a GH(16,1)
//D3. Take each row of GH(16,1), repeat it 16 times and put it in a matrix.In this way, 16 matrices are obtained. The ordered collection of these  matrices is called R[s][s][s].
//D4. Find the Hadamard product of each of the above matrices with a Hadamard matrix of order 16.
//D5. Thus 16 matrices of the order 16x16 have been obtained. 
//D6. Multiply each of the above matrices with a diagonal matrix consisting of a Hadamard matrix H.
//7. Similarly multiply each matrix with a diagonal matrix that contains one of many in a class of unbiased Hadamard matrices of order 16.
//D8. Convert each member of the group to a row transformation matrix applied on the Hadamard matrix of order 16.
//D9. This should give you another Hadamard matrix.
//8. Count the number of and check if the resulting Hadamard matrices are unbiased.

//Global variable
//integers
const int s=16;//size of our matrix in the form of 2^k
const int s2 = s*s;
const int k=4;//k
const int number_of_unbiased_Hadamard_of_order_s=s/2;
//matrix

int unbiased_Hadamard_of_order_s[s/2][s][s]={{{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},{1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1},{1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1},{1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1},{1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1},{1,-1,1,-1,-1,1,-1,1,1,-1,1,-1,-1,1,-1,1},{1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1},{1,-1,-1,1,-1,1,1,-1,1,-1,-1,1,-1,1,1,-1},{1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1},{1,-1,1,-1,1,-1,1,-1,-1,1,-1,1,-1,1,-1,1},{1,1,-1,-1,1,1,-1,-1,-1,-1,1,1,-1,-1,1,1},{1,-1,-1,1,1,-1,-1,1,-1,1,1,-1,-1,1,1,-1},{1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1},{1,-1,1,-1,-1,1,-1,1,-1,1,-1,1,1,-1,1,-1},{1,1,-1,-1,-1,-1,1,1,-1,-1,1,1,1,1,-1,-1},{1,-1,-1,1,-1,1,1,-1,-1,1,1,-1,1,-1,-1,1}},{{1,1,1,1,1,-1,1,-1,1,1,-1,-1,-1,1,1,-1},{1,-1,1,-1,1,1,1,1,1,-1,-1,1,-1,-1,1,1},{1,1,-1,-1,1,-1,-1,1,1,1,1,1,-1,1,-1,1},{1,-1,-1,1,1,1,-1,-1,1,-1,1,-1,-1,-1,-1,-1},{1,1,1,1,1,-1,1,-1,-1,-1,1,1,1,-1,-1,1},{1,-1,1,-1,1,1,1,1,-1,1,1,-1,1,1,-1,-1},{1,1,-1,-1,1,-1,-1,1,-1,-1,-1,-1,1,-1,1,-1},{1,-1,-1,1,1,1,-1,-1,-1,1,-1,1,1,1,1,1},{1,1,1,1,-1,1,-1,1,1,1,-1,-1,1,-1,-1,1},{1,-1,1,-1,-1,-1,-1,-1,1,-1,-1,1,1,1,-1,-1},{1,1,-1,-1,-1,1,1,-1,1,1,1,1,1,-1,1,-1},{1,-1,-1,1,-1,-1,1,1,1,-1,1,-1,1,1,1,1},{-1,-1,-1,-1,1,-1,1,-1,1,1,-1,-1,1,-1,-1,1},{-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,1,-1,-1},{-1,-1,1,1,1,-1,-1,1,1,1,1,1,1,-1,1,-1},{-1,1,1,-1,1,1,-1,-1,1,-1,1,-1,1,1,1,1}},{{1,1,1,1,1,1,-1,-1,1,-1,-1,1,-1,1,-1,1},{1,-1,1,-1,1,-1,-1,1,1,1,-1,-1,-1,-1,-1,-1},{1,1,-1,-1,1,1,1,1,1,-1,1,-1,-1,1,1,-1},{1,-1,-1,1,1,-1,1,-1,1,1,1,1,-1,-1,1,1},{1,1,1,1,1,1,-1,-1,-1,1,1,-1,1,-1,1,-1},{1,-1,1,-1,1,-1,-1,1,-1,-1,1,1,1,1,1,1},{1,1,-1,-1,1,1,1,1,-1,1,-1,1,1,-1,-1,1},{1,-1,-1,1,1,-1,1,-1,-1,-1,-1,-1,1,1,-1,-1},{1,1,1,1,-1,-1,1,1,1,-1,-1,1,1,-1,1,-1},{1,-1,1,-1,-1,1,1,-1,1,1,-1,-1,1,1,1,1},{1,1,-1,-1,-1,-1,-1,-1,1,-1,1,-1,1,-1,-1,1},{1,-1,-1,1,-1,1,-1,1,1,1,1,1,1,1,-1,-1},{-1,-1,-1,-1,1,1,-1,-1,1,-1,-1,1,1,-1,1,-1},{-1,1,-1,1,1,-1,-1,1,1,1,-1,-1,1,1,1,1},{-1,-1,1,1,1,1,1,1,1,-1,1,-1,1,-1,-1,1},{-1,1,1,-1,1,-1,1,-1,1,1,1,1,1,1,-1,-1}},{{1,1,1,1,1,-1,-1,1,1,-1,1,-1,-1,-1,1,1},{1,-1,1,-1,1,1,-1,-1,1,1,1,1,-1,1,1,-1},{1,1,-1,-1,1,-1,1,-1,1,-1,-1,1,-1,-1,-1,-1},{1,-1,-1,1,1,1,1,1,1,1,-1,-1,-1,1,-1,1},{1,1,1,1,1,-1,-1,1,-1,1,-1,1,1,1,-1,-1},{1,-1,1,-1,1,1,-1,-1,-1,-1,-1,-1,1,-1,-1,1},{1,1,-1,-1,1,-1,1,-1,-1,1,1,-1,1,1,1,1},{1,-1,-1,1,1,1,1,1,-1,-1,1,1,1,-1,1,-1},{1,1,1,1,-1,1,1,-1,1,-1,1,-1,1,1,-1,-1},{1,-1,1,-1,-1,-1,1,1,1,1,1,1,1,-1,-1,1},{1,1,-1,-1,-1,1,-1,1,1,-1,-1,1,1,1,1,1},{1,-1,-1,1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,-1},{-1,-1,-1,-1,1,-1,-1,1,1,-1,1,-1,1,1,-1,-1},{-1,1,-1,1,1,1,-1,-1,1,1,1,1,1,-1,-1,1},{-1,-1,1,1,1,-1,1,-1,1,-1,-1,1,1,1,1,1},{-1,1,1,-1,1,1,1,1,1,1,-1,-1,1,-1,1,-1}},
	{{1,1,1,-1,1,1,1,-1,1,1,1,-1,-1,-1,-1,1},{1,1,-1,1,1,1,-1,1,1,1,-1,1,-1,-1,1,-1},{1,-1,1,1,1,-1,1,1,1,-1,1,1,-1,1,-1,-1},{-1,1,1,1,-1,1,1,1,-1,1,1,1,1,-1,-1,-1},{1,1,1,-1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1},{1,1,-1,1,1,1,-1,1,-1,-1,1,-1,1,1,-1,1},{1,-1,1,1,1,-1,1,1,-1,1,-1,-1,1,-1,1,1},{-1,1,1,1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1},{1,1,1,-1,-1,-1,-1,1,1,1,1,-1,1,1,1,-1},{1,1,-1,1,-1,-1,1,-1,1,1,-1,1,1,1,-1,1},{1,-1,1,1,-1,1,-1,-1,1,-1,1,1,1,-1,1,1},{-1,1,1,1,1,-1,-1,-1,-1,1,1,1,-1,1,1,1},{-1,-1,-1,1,1,1,1,-1,1,1,1,-1,1,1,1,-1},{-1,-1,1,-1,1,1,-1,1,1,1,-1,1,1,1,-1,1},{-1,1,-1,-1,1,-1,1,1,1,-1,1,1,1,-1,1,1},{1,-1,-1,-1,-1,1,1,1,-1,1,1,1,-1,1,1,1}},{{1,1,1,-1,1,1,-1,1,1,-1,1,1,1,-1,-1,-1},{1,1,-1,1,1,1,1,-1,-1,1,1,1,-1,1,-1,-1},{1,-1,1,1,-1,1,1,1,1,1,1,-1,-1,-1,1,-1},{-1,1,1,1,1,-1,1,1,1,1,-1,1,-1,-1,-1,1},{1,1,1,-1,1,1,-1,1,-1,1,-1,-1,-1,1,1,1},{1,1,-1,1,1,1,1,-1,1,-1,-1,-1,1,-1,1,1},{1,-1,1,1,-1,1,1,1,-1,-1,-1,1,1,1,-1,1},{-1,1,1,1,1,-1,1,1,-1,-1,1,-1,1,1,1,-1},{1,1,1,-1,-1,-1,1,-1,1,-1,1,1,-1,1,1,1},{1,1,-1,1,-1,-1,-1,1,-1,1,1,1,1,-1,1,1},{1,-1,1,1,1,-1,-1,-1,1,1,1,-1,1,1,-1,1},{-1,1,1,1,-1,1,-1,-1,1,1,-1,1,1,1,1,-1},{-1,-1,-1,1,1,1,-1,1,1,-1,1,1,-1,1,1,1},{-1,-1,1,-1,1,1,1,-1,-1,1,1,1,1,-1,1,1},{-1,1,-1,-1,-1,1,1,1,1,1,1,-1,1,1,-1,1},{1,-1,-1,-1,1,-1,1,1,1,1,-1,1,1,1,1,-1}},{{1,1,1,-1,1,-1,1,1,-1,1,1,1,-1,-1,1,-1},{1,1,-1,1,-1,1,1,1,1,-1,1,1,-1,-1,-1,1},{1,-1,1,1,1,1,1,-1,1,1,-1,1,1,-1,-1,-1},{-1,1,1,1,1,1,-1,1,1,1,1,-1,-1,1,-1,-1},{1,1,1,-1,1,-1,1,1,1,-1,-1,-1,1,1,-1,1},{1,1,-1,1,-1,1,1,1,-1,1,-1,-1,1,1,1,-1},{1,-1,1,1,1,1,1,-1,-1,-1,1,-1,-1,1,1,1},{-1,1,1,1,1,1,-1,1,-1,-1,-1,1,1,-1,1,1},{1,1,1,-1,-1,1,-1,-1,-1,1,1,1,1,1,-1,1},{1,1,-1,1,1,-1,-1,-1,1,-1,1,1,1,1,1,-1},{1,-1,1,1,-1,-1,-1,1,1,1,-1,1,-1,1,1,1},{-1,1,1,1,-1,-1,1,-1,1,1,1,-1,1,-1,1,1},{-1,-1,-1,1,1,-1,1,1,-1,1,1,1,1,1,-1,1},{-1,-1,1,-1,-1,1,1,1,1,-1,1,1,1,1,1,-1},{-1,1,-1,-1,1,1,1,-1,1,1,-1,1,-1,1,1,1},{1,-1,-1,-1,1,1,-1,1,1,1,1,-1,1,-1,1,1}},{{1,1,1,-1,-1,1,1,1,1,1,-1,1,-1,1,-1,-1},{1,1,-1,1,1,-1,1,1,1,1,1,-1,1,-1,-1,-1},{1,-1,1,1,1,1,-1,1,-1,1,1,1,-1,-1,-1,1},{-1,1,1,1,1,1,1,-1,1,-1,1,1,-1,-1,1,-1},{1,1,1,-1,-1,1,1,1,-1,-1,1,-1,1,-1,1,1},{1,1,-1,1,1,-1,1,1,-1,-1,-1,1,-1,1,1,1},{1,-1,1,1,1,1,-1,1,1,-1,-1,-1,1,1,1,-1},{-1,1,1,1,1,1,1,-1,-1,1,-1,-1,1,1,-1,1},{1,1,1,-1,1,-1,-1,-1,1,1,-1,1,1,-1,1,1},{1,1,-1,1,-1,1,-1,-1,1,1,1,-1,-1,1,1,1},{1,-1,1,1,-1,-1,1,-1,-1,1,1,1,1,1,1,-1},{-1,1,1,1,-1,-1,-1,1,1,-1,1,1,1,1,-1,1},{-1,-1,-1,1,-1,1,1,1,1,1,-1,1,1,-1,1,1},{-1,-1,1,-1,1,-1,1,1,1,1,1,-1,-1,1,1,1},{-1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1},{1,-1,-1,-1,1,1,1,-1,1,-1,1,1,1,1,-1,1}}};
int inner_product(int* a, int*b, int n);
int** make_Hadamard_matrix(int order);
int** negative_matrix(int** a, int n);
//This function generates a Hadamard matrix by multiplying a diagonal row of Hadamard matrices unbiased_Hadamard_of_order_s[n] with the matrix M[m]
int** unbiased_Hadamard_of_order_s2(int n, int m); 
bool check_Hadamard(int**, int);

int product[s][s];//multiplication table will save in it if we call the function multiplication_table
int H[s][s];//this is a hadamard matrix of size s in the form we use in our report. 
int R[s][s][s];//Made by reapition of rows of generalized Hadamard matrix
int IR[s][s][s]={0};// this is the matrix of form IxIxI, IxIxR, .... we should call the function generate_IR to fill this function bace in index i
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
void group_operation_table();//save the multiplication table in matrix product using multiplication function
int group_operation(int ,int);//return the multiplication of two arguments using binary code
int*  binary(int );//change the argument to binary
int integer(int* );//recieve a binary and change it to an integer

//the six following functions do the normal matrix multiplication 
int** matrix_multiplication(int[s][s] ,int [s][s]);
int** matrix_multiplication(int[s][s][s] ,int [s][s][s]);
int** matrix_multiplication(int**,int [s][s]);
int** matrix_multiplication(int[s][s] ,int**);
int** matrix_multiplication(int a[s][s2][s2],int** b, int y);
int** matrix_multiplication(int** ,int**, int);
////return the transposition of the argument
int** Trans(int a[s][s]);
int** Trans(int a[s][s2][s2], int);
int** Trans(int**a , int n);
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
void generate_IR (int,int ,int ,int ,int );//save the IxIxR in IR matrix, the third int shoul be 1.


void fill_row();//this function fill the variable row
void generate_GH();//this function save a generalized hadamard matrix in variable p[][] which is anarray object of class polynomail
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


//the two following variables are used to make a generalized hadamard matrix
polynomial p[s][s];
polynomial row[s];

int main()
{ 
   
  fill_row();//this function fill the variable row
  group_operation_table();//save the multiplication table in matrix product using multiplication function
  generate_GH();//this function save a generalized hadamard matrix in variable p[][] which is anarray object of class polynomail
  for(int i=0; i < s; i++)
	 {      
		 IR[i][0][i]=1;
       	         generate_IR(0,i,1,i,i);
	 }
			        
  make_matrix_R();//this function fill R[s][s][s]
  make_matrix_M();//fill w by componentwise product of H and  R
 // print(unbiased_Hadamard_of_order_s2(1,1),s2);
 bool check[number_of_unbiased_Hadamard_of_order_s*s][number_of_unbiased_Hadamard_of_order_s*s]={{false}};
bool unbiased[number_of_unbiased_Hadamard_of_order_s*s][number_of_unbiased_Hadamard_of_order_s*s]={{false}};
 int index_i=-1,index_j=-1;

  for(int n=0;n<number_of_unbiased_Hadamard_of_order_s;n++)
  {
	  for(int m=0;m<s;m++)
	  {
		  index_i++;
		  index_j=-1;
		  for(int i=0;i<number_of_unbiased_Hadamard_of_order_s;i++)
		  {
			  
			  for(int j=0;j<s;j++)
			  {
				index_j++;
				if(check[index_i][index_j]!= true && index_i!=index_j)
				  {
					  
					  unbiased[index_i][index_j]=check_Hadamard(matrix_multiplication(unbiased_Hadamard_of_order_s2(n,m),Trans(unbiased_Hadamard_of_order_s2(i,j),s2),s2),s2);
					  check[index_i][index_j]=true;
					  check[index_j][index_i]=true;
					  unbiased[index_j][index_i] = unbiased[index_i][index_j];
				  }
				if(index_i==index_j)
					unbiased[index_i][index_j]=true;

				  cout<<unbiased[index_i][index_j];
			  }
		  }
	cout<<endl;	  
	  }
	  
  }

      
/*   for(int i=0;i<s;i++)
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
*/
 int a=1,b=1,c=1,d=2; 
cout<<check_Hadamard(matrix_multiplication(unbiased_Hadamard_of_order_s2(a,b),Trans(unbiased_Hadamard_of_order_s2(c,d),s2),s2),s2);
  return 0;

}

/*
bool check_if_2_matrices_are_equal(int** a, int** b,int n)
{
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			if(a[i][j]!=b[i][j])
				return false;
		}
	}
	return true;
}
*/
int** Trans(int** a, int n)
{
	int** trans = new int*[n];
	for(int i=0;i<n;i++)
	{
		trans[i] = new int[n];
	}
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			trans[i][j]=a[j][i];
		}
	}
	return trans;
}

bool check_Hadamard(int** a, int n)
{
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			if(a[i][j]!=16 && a[i][j]!=-16)
				return false;
		}
	}
	for(int i=0;i<n;i++)
	{
		for(int j=i+1;j<n;j++)
		{
			if(inner_product(a[i],a[j],n)!=0)
			{
				return false;
			}
		}
	}
	return true;
}

//Generates a Hadamard matrix of order order(has to be a power of 2). 
int** make_Hadamard_matrix(int order)
{
	int** a = new int*[s];
        for(int i=0;i<s;i++)
	{
		a[i] = new int[s];
	}
	for(int i=0;i<s;i++)
	{
		for(int j=0;j<s;j++)
		{
			a[i][j] = pow(-1,inner_product(binary(i),binary(j),k));
		}
	}
	
	return (a);

}


//Returns inner product of 2 arrays.
int inner_product(int* a,int* b, int n)
{
  int sum=0; 
  for(int i=0;i<n;i++)
   {
	sum = sum + a[i]*b[i];
   }
  return(sum);
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

//Generates all rows of length index, and checks if any of them are unbiased with the row_index^th row of the a_index^th matrix of the 3-dimensional matrix a. Returns true if the row is unbiased 
//and false otherwise.
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

// This function prints the negative of a number.
int** negative_matrix(int** a, int n)
{
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			a[i][j] = -a[i][j];
		}
	}
	return a;
}

//This function generates a Hadamard matrix by multiplying a diagonal row of Hadamard matrices unbiased_Hadamard_of_order_s[n] with the matrix M[m]
int** unbiased_Hadamard_of_order_s2(int n, int m)
{
	
	int** Hadamard = new int*[s2];
	for(int i=0;i<s2;i++)
	{
		Hadamard[i] = new int[s2];
	}
	int** p ;
	for(int i=0;i<s;i++)
	{
		for(int j=0;j<s;j++)
		{
			//print(matrix_multiplication(unbiased_Hadamard_of_order_s[n], IR[abs(M[m][i][j])]));
		//	print(IR[abs(M[m][i][j])]);
			//print(unbiased_Hadamard_of_order_s[n]);
			
			p = matrix_multiplication(unbiased_Hadamard_of_order_s[n], IR[abs(M[m][i][j])-1]);
			if(M[m][i][j]<0)
			{
				p = negative_matrix(p,s);
			}
			
			//cout<<endl;
			//cout<<"Hadamard matrix ="<<endl;
			//print(unbiased_Hadamard_of_order_s[n]);
			//cout<<"IR matrix ="<<endl;
			//print(IR[abs(M[m][i][j])-1]);
			//cout<<"Product matrix p="<<endl;
			//print(p);	
			for(int k=0;k<s;k++)
			{
				for(int l=0;l<s;l++)
				{
					Hadamard[i*s +k][j*s +l]=p[k][l];
				}
			}
		}
	}
	return Hadamard;
}

//This function only checks if the matrix expanded from M[d], which is in replac[d] is a weighing matrix or not. 
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

// Checks if a matrix A is of Bush type or not.
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

//Computes the operation table of the group G, over which the generalized Hadamard matrix is constructed. The elements of G are indexed by numbers, each element of G is the binary expansion of 
//that index.
void group_operation_table()
{
   for(int i= 0; i < s; i++)
   {
      for(int j = 0; j < s; j++)
      {
	 product[i][j] = group_operation(i,j);
      }
   }
 }

//The group operation of G over which the Generalized Hadamard matrix is construct has as its group operation, the addition of binary representations of the 2 elements. The parameters of the
//function are indices of the group elements whose operation is to be computed. 
int group_operation(int y, int z)
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

//Returns the binary array of an integer a.
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

//Returns the integer corresponding to the binary expansion in the array a.
int integer(int* a)
{
  
   int result=0;
   for(int i= 0; i< k; i++)
     {
	result = result + (a[i]*pow(2.0,i));

     }
    return result;
}

//The following matrix_multiplication functions multiply 2 matrices in the usual sense.
int** matrix_multiplication(int a[s][s],int b[s][s])
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
int** matrix_multiplication(int a[s][s2][s2] ,int b[s][s2][s2], int y)
{
	cout<<"2"<<endl;
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
   
	cout<<"3"<<endl;
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
   
	cout<<"4"<<endl;
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
int** matrix_multiplication(int** a,int** b,int n)
{
   int** c =new int*[n];
   for(int i=0; i < n; i++)
   {
      c[i] = new int[n];
   }
   int sum;
   for(int i = 0; i < n; i++)
   {
      for(int j = 0; j < n; j++)
      {
	 sum =0;
	 for(int o =0; o < n; o++)
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

void print(int b[s2][s2])
{
	for(int i=0;i<s2;i++)
	{
		for(int j=0;j<s2;j++)
		{
			cout<<b[i][j]<<" ";
		}
		cout<<endl;
	}
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

void generate_IR(int indexi,int indexj,int size,int i,int index_IR)
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
		generate_IR(indexi,indexj,2*size,i,index_IR);
        }
        else if(q==1)
        {
                for(int m=0;m<size;m++)
                        for(int l=0;l<size;l++)
                        {
			   IR[index_IR][indexi+size+m][indexj-size+l]=IR[index_IR][indexi+m][indexj+l];
		
                        }
		generate_IR(indexi,indexj-size,2*size,i,index_IR);
        }
}

/*
 * I dont think this function is needed.
void generate_HH()
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
*/

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




void generate_GH()
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
 int** H;
	H=make_Hadamard_matrix(s);
for(int i=0;i<s;i++)
{
	for(int j=0;j<s;j++)
	{
		for(int l=0;l<s;l++)
		{
			M[i][j][l] = H[j][l]*R[i][j][l];
		}
	}
}
/*	
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
*/
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
   print(R);
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
