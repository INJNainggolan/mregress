#include<iostream>
#include<fstream>
#include<iomanip>
using namespace std;

void transpose(double **p1,double **p2,int m,int n);
void multipl(double **p1,double **p2,double **p3,int m,int n,int p);
void Inver(double **p1,double **p2,int n);
double SD(double **p1,double **p2,double **p3,double **p4,int m,int n);
double ST(double **p1,int m);
void de_allocate(double **data,int m);


int main() {
int row,col;
char filename[30];
double SDsum,STsum,F,R2;

cout<<"Input original data file: \n";
ifstream infile;  //���ļ�
cin>>filename;
infile.open(filename);
if(!infile) {
cout<<"Opening the file failed!\n";
exit(1);
}
infile>>row>>col; //�����ļ��е�����������

double **matrix=new double*[row]; //Ϊ��̬��ά��������ڴ�
double **X=new double*[row];
double **Y=new double*[row];
double **XT=new double*[col];
double **XTX=new double*[col];
double **XTXInv=new double*[col];
double **XTXInvXT=new double*[col];
double **B=new double*[col];
double **YE=new double*[row];
for(int i=0;i<row;i++) {
  matrix[i]=new double[col];
  X[i]=new double[col];
  Y[i]=new double[1];
  Y[i]=new double[1];
  YE[i]=new double[1];
}
for(int i=0;i<col;i++) {
  XT[i]=new double[row];
  XTX[i]=new double[2*col];/////////////////////Ϊʲô�������2*col�пռ������col���ھ�������ʱ��XTX���������������Ϊԭ��2�ɱ����������㷨�йء�
  XTXInv[i]=new double[col];
  XTXInvXT[i]=new double[row];
  B[i]=new double[1];
}
for(int i=0;i<row;i++)
  for(int j=0;j<col;j++)
    infile>>matrix[i][j];
infile.close();

for(int i=0;i<row;i++) { //��ȡ1X��Y������
  X[i][0]=1;
  Y[i][0]=matrix[i][col-1];
  for(int j=0;j<col-1;j++)
    X[i][j+1]=matrix[i][j];
}

transpose(X,XT,row,col);
multipl(XT,X,XTX,col,row,col);
Inver(XTX,XTXInv,col);
multipl(XTXInv,XT,XTXInvXT,col,col,row);
multipl(XTXInvXT,Y,B,col,row,1);
SDsum=SD(Y,X,B,YE,row,col);
STsum=ST(Y,row);
F=((STsum-SDsum)/(col-1))/(SDsum/(row-col));
R2=1/(1+(row-col)/F/(col-1));


cout<<"���B:\n";  //��Ļ������B��SD��ST��F��R2
for(int i=0;i<col;i++)
  cout<<setiosflags(ios::fixed)<<setprecision(4)<<B[i][0]<<' ';
  cout<<endl;
cout<<"SD="<<SDsum<<';'<<"ST="<<STsum<<';'<<"F="<<F<<';'<<"R2="<<R2<<endl;


ofstream outfile; // ���д���ļ�
cout<<"Output file'name:\n";
cin>>filename;
outfile.open(filename);
if(!outfile) {
  cout<<"Opening the file failed!\n";
  exit(1);
}
outfile<<"���B:\n";
for(int i=0;i<col;i++)
  outfile<<B[i][0]<<' ';
  outfile<<endl;
outfile<<setiosflags(ios::fixed)<<setprecision(4)<<"SD="<<SDsum<<';'<<"ST="<<STsum<<';'<<"F="<<F<<';'<<"R2="<<R2<<endl;
outfile<<"Y and YE and Y-YE's value are:\n";
for(int i=0;i<row;i++)
	outfile<<Y[i][0]<<"  "<<YE[i][0]<<"  "<<Y[i][0]-YE[i][0]<<endl;
outfile.close();


de_allocate(matrix,row);
de_allocate(X,row);
de_allocate(Y,row);
de_allocate(XT,col);
de_allocate(XTX,col);
de_allocate(XTXInv,col);
de_allocate(XTXInvXT,col);
de_allocate(B,col);
de_allocate(YE,row);

system("pause");
return(0);
}

void de_allocate(double **data,int m) { //�ͷ��ڴ浥Ԫ
  for(int i=0;i<m;i++)
    delete []data[i];
  delete []data;
}

double ST(double **p1,int m) { //�������ƽ����ST
  double sum1=0,sum2=0,Yave=0;
  for(int i=0;i<m;i++)
    sum1+=p1[i][0];
  Yave=sum1/m;
  for(int i=0;i<m;i++)
    sum2+=(p1[i][0]-Yave)*(p1[i][0]-Yave);
  return sum2;
}

double SD(double **p1,double **p2,double **p3,double **p4,int m,int n) { //��ƫ��ƽ����SD
  double sum1=0,sum2=0;
  for(int i=0;i<m;i++) {
	sum1=0;
    for(int k=0;k<n;k++)
	  sum1+=p2[i][k]*p3[k][0];
    p4[i][0]=sum1;
  }
  for(int i=0;i<m;i++)
    sum2+=(p1[i][0]-p4[i][0])*(p1[i][0]-p4[i][0]);
  return sum2;
}


void transpose(double **p1,double **p2,int m,int n) {  //����ת��
for(int i=0;i<n;i++)
  for(int j=0;j<m;j++)
      p2[i][j]=p1[j][i];
}

void multipl(double **p1,double **p2,double **p3,int m,int n,int p) { //�������
double sum;
for(int i=0;i<m;i++) {
  for(int j=0;j<p;j++) {
	sum=0;
    for(int k=0;k<n;k++)
	  sum+=p1[i][k]*p2[k][j];
    p3[i][j]=sum;
  }
}
}

void Inver(double **p1,double **p2,int n) {  //�������

//��ʼ���������Ҳ���뵥λ��
  for(int i=0;i<n;i++) {
	for(int j=0;j<n;j++) {
	  p1[i][j+n]=0;
	  p1[i][i+n]=1;
	}
  }

//���ڶԽ�Ԫ��Ϊ0�Ľ��л��в���
  for(int i=0;i<n;i++) 
   {
    while(p1[i][i]==0)
	 {
      for(int j=i+1;j<n;j++)
	    {
	     if (p1[j][i]!=0) 
		 { double temp=0;
		  for(int r=i;r<2*n;r++)
		    {temp=p1[j][r];p1[j][r]=p1[i][r];p1[i][r]=temp;}
		 }
	     break;
	    }
	 }
    //if (p1[i][i]==0) return 0;
   }
 //�б任Ϊ�����Ǿ���
  double k=0;
  for(int i=0;i<n;i++) {
	for(int j=i+1;j<n;j++) {
	k=(-1)*p1[j][i]/p1[i][i];
    for(int r=i;r<2*n;r++)
	  p1[j][r]+=k*p1[i][r];
	}
  }
  //�б任Ϊ�����Ǿ���
  //double k=0;
  for(int i=n-1;i>=0;i--) {
	for(int j=i-1;j>=0;j--) {
	  k=(-1)*p1[j][i]/p1[i][i];
	  for(int r=0;r<2*n;r++)
        p1[j][r]+=k*p1[i][r];
	  }
  }
  //��Ϊ��λ��
  for(int i=n-1;i>=0;i--) {
    k=p1[i][i];
	for(int j=0;j<2*n;j++) 
	  p1[i][j]/=k;
  }

  //��ֳ������
  for(int i=0;i<n;i++) {
    for(int j=0;j<n;j++)
	  p2[i][j]=p1[i][n+j];
  }

}


