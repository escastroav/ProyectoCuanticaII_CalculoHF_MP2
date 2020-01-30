#include <cmath>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

//Definimos las variables y matrices que todo el programa va a usar
Matrix2d T,NucA,NucB,H,S,vecS,autS,Xcanon,G,C,P,Oldp,F,Fp,Cp,E;
VectorXd HFx_axis(81),HFener(81),MP2x_axis(81),MP2corr(81),MP2ener(81);
double Delta = 0;
double Energy = 0;
//Definicion de funciones
void Init2d(double TT[2][2][2][2])
{
  for(int i=0;i<2;i++){
    for(int j=0;j<2;j++){
      for(int k=0;k<2;k++){
	for(int l=0;l<2;l++){
	  TT[i][j][k][l] = 0;
	}
      }
    }
  }
}
void Init4d(double TT[4][4][4][4])
{
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      for(int k=0;k<4;k++){
	for(int l=0;l<4;l++){
	  TT[i][j][k][l] = 0;
	}
      }
    }
  }
}
double TT[2][2][2][2];
double S_int(double A, double B, double Rab2)
{
  return pow(M_PI/(A+B),1.5)*exp(-A*B*Rab2/(A+B));
}
double T_int(double A, double B, double Rab2)
{
  return A*B/(A+B)*(3.0-2.0*A*B*Rab2/(A+B))*pow(M_PI/(A+B),1.5)*exp(-A*B*Rab2/(A+B));
}
double F0(double t)
{
  if(t<1e-6)
    return 1.0-t/3.0;
  else
    return 0.5*pow(M_PI/t,0.5)*erf(sqrt(t));
}
double V_int(double A, double B, double Rab2, double Rcp2, double Zc)
{
  double V = 2.0*M_PI/(A+B)*F0((A+B)*Rcp2)*exp(-A*B*Rab2/(A+B));
  return -V*Zc;
}
double TwoE(double A, double B, double C, double D, double Rab2, double Rcd2, double Rpq2)
{
  return 2.0*pow(M_PI,2.5)/((A+B)*(C+D)*sqrt(A+B+C+D))*F0((A+B)*(C+D)*Rpq2/(A+B+C+D))*exp(-A*B*Rab2/(A+B)-C*D*Rcd2/(C+D));
}
double S12,T11,T12,T22,V11A,V12A,V22A,V11B,V12B,V22B,V1111,V2111,V2121,V2211,V2221,V2222;
double Rap = 0,Rbp = 0,Raq = 0,Rbq = 0,Rpq = 0,Rap2 = 0,Rbp2 = 0,Raq2 = 0,Rbq2 = 0,Rpq2 = 0;
void Intgrl(int N,double R,double Zeta1, double Zeta2, double Za, double Zb)
{
  S12 = 0;
  T11 = 0;
  T12 = 0;
  T22 = 0;
  V11A = 0;
  V12A = 0;
  V22A = 0;
  V11B = 0;
  V12B = 0;
  V22B = 0;
  V1111 = 0;
  V2111 = 0;
  V2121 = 0;
  V2211 = 0;
  V2221 = 0;
  V2222 = 0;
  
    
  double R2 = R*R;

  Matrix3d Coeff;  Coeff << 1.000000,0.000000,0.000000,
			0.678914,0.430129,0.000000,
			0.444635,0.535328,0.154329;

  Matrix3d Expon; Expon << 0.270950,0.000000,0.000000,
			0.151623,0.851819,0.000000,
			0.109818,0.405771,2.227660;
  VectorXd D1(3),
    A1(3),
    D2(3),
    A2(3);

  for(int i=0;i<N;i++){
    A1(i) = Expon(N-1,i)*pow(Zeta1,2); //Exponente
    D1(i) = Coeff(N-1,i)*pow(2.0*A1(i)/M_PI,0.75); //Norm
    A2(i) = Expon(N-1,i)*pow(Zeta2,2);
    D2(i) = Coeff(N-1,i)*pow(2.0*A2(i)/M_PI,0.75);  
  }
  
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      Rap = A2(j)*R/(A1(i)+A2(j));
      Rap2 = pow(Rap,2);
      Rbp2 = pow(R-Rap,2);
      S12 = S12 + S_int(A1(i),A2(j),R2)*D1(i)*D2(j);
      T11 = T11 + T_int(A1(i),A1(j),0.0)*D1(i)*D1(j);
      T12 = T12 + T_int(A1(i),A2(j),R2)*D1(i)*D2(j);
      T22 = T22 + T_int(A2(i),A2(j),0.0)*D2(i)*D2(j);
      V11A = V11A + V_int(A1(i),A1(j),0.0,0.0,Za)*D1(i)*D1(j);
      V12A = V12A + V_int(A1(i),A2(j),R2,Rap2,Za)*D1(i)*D2(j);
      V22A = V22A + V_int(A2(i),A2(j),0.0,R2,Za)*D2(i)*D2(j);
      V11B = V11B + V_int(A1(i),A1(j),0.0,R2,Zb)*D1(i)*D1(j);
      V12B = V12B + V_int(A1(i),A2(j),R2,Rbp2,Zb)*D1(i)*D2(j);
      V22B = V22B + V_int(A2(i),A2(j),0.0,0.0,Zb)*D2(i)*D2(j);
    }
  }
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      for(int k=0;k<N;k++){
	for(int l=0;l<N;l++){
	  Rap = A2(i)*R/(A2(i)+A1(j));
	  Rbp = R - Rap;
	  Raq = A2(k)*R/(A2(k)+A1(l));
	  Rbq = R - Raq;
	  Rpq = Rap - Raq;
	  Rap2 = Rap*Rap;
	  Rbp2 = Rbp*Rbp;
	  Raq2 = Raq*Raq;
	  Rbq2 = Rbq*Rbq;
	  Rpq2 = Rpq*Rpq;
	  V1111 = V1111 + TwoE(A1(i),A1(j),A1(k),A1(l),0.0,0.0,0.0)*D1(i)*D1(j)*D1(k)*D1(l);
	  V2111 = V2111 + TwoE(A2(i),A1(j),A1(k),A1(l),R2,0.0,Rap2)*D2(i)*D1(j)*D1(k)*D1(l);
	  V2121 = V2121 + TwoE(A2(i),A1(j),A2(k),A1(l),R2,R2,Rpq2)*D2(i)*D1(j)*D2(k)*D1(l);
	  V2211 = V2211 + TwoE(A2(i),A2(j),A1(k),A1(l),0.0,0.0,R2)*D2(i)*D2(j)*D1(k)*D1(l);
	  V2221 = V2221 + TwoE(A2(i),A2(j),A2(k),A1(l),0.0,R2,Rbq2)*D2(i)*D2(j)*D2(k)*D1(l);
	  V2222 = V2222 + TwoE(A2(i),A2(j),A2(k),A2(l),0.0,0.0,0.0)*D2(i)*D2(j)*D2(k)*D2(l);
	}
      }
    }
  }
  //cout << V1111 << endl;
  //cout << V2111 << endl;
  //cout << V2121 << endl;
  //cout << V2211 << endl;
  //cout << V2221 << endl;
  //cout << V2222 << endl;
}

void DiagS(Matrix2d A)
{
  SelfAdjointEigenSolver<Matrix2d> eigenSolver(A);
  Vector2d values = eigenSolver.eigenvalues();
  vecS = eigenSolver.eigenvectors();
  autS = values.asDiagonal();
}

void DiagFp(Matrix2d A)
{
  SelfAdjointEigenSolver<Matrix2d> eigenSolver(A);
  Vector2d values = eigenSolver.eigenvalues();
  Cp = eigenSolver.eigenvectors();
  E = values.asDiagonal();
}

void Colect(int N, double R, double Zeta1, double Zeta2, double Za, double Zb)
{
 //Matriz Cinetica
  T(0,0) = T11;
  T(0,1) = T12;
  T(1,0) = T12;
  T(1,1) = T22;
  
  //Matriz nuclear A
  NucA(0,0) = V11A;
  NucA(0,1) = V12A;
  NucA(1,0) = V12A;
  NucA(1,1) = V22A;
  
  //Matriz nuclear B
  NucB(0,0) = V11B;
  NucB(0,1) = V12B;
  NucB(1,0) = V12B;
  NucB(1,1) = V22B;
  
  // Matriz de HCore
  H(0,0) = T11+V11A+V11B;
  H(0,1) = T12+V12A+V12B;
  H(1,0) = H(0,1);
  H(1,1) = T22+V22A+V22B;
  
  // Matriz de Overlap
  S(0,0) = 1.0;
  S(0,1) = S12;
  S(1,0) = S12;
  S(1,1) = 1.0;
  
  //Definicion de X
  DiagS(S);
  Matrix2d smum = autS.array().sqrt().inverse().matrix().diagonal().asDiagonal();
  Matrix2d Smum = vecS * (smum * vecS.transpose());
  Xcanon = vecS * smum;

  //Matriz coulombiana y de intercambio
  TT[0][0][0][0] = V1111;
  TT[1][0][0][0] = V2111;
  TT[0][1][0][0] = V2111;
  TT[0][0][1][0] = V2111;
  TT[0][0][0][1] = V2111;
  TT[1][0][1][0] = V2121;
  TT[0][1][1][0] = V2121;
  TT[1][0][0][1] = V2121;
  TT[0][1][0][1] = V2121;
  TT[1][1][0][0] = V2211;
  TT[0][0][1][1] = V2211;
  TT[1][1][1][0] = V2221;
  TT[1][1][0][1] = V2221;
  TT[1][0][1][1] = V2221;
  TT[0][1][1][1] = V2221;
  TT[1][1][1][1] = V2222;
}
double Energytot = 0;
void SCF(int N, double R, double Zeta1, double Zeta2, double Za, double Zb, Matrix2d G)
{
  //cout << Xcanon << endl;
  double Crit = 1e-11;//Criterio de convergencia
  int Maxit = 250;//Numero maximo de iteraciones
  int Iter = 0;

  //----- Paso 1. El guess inicial corresponde a P=0 -----//
  // El primer calculo se hara sin repulsiones electronicas
  
  P = Matrix2d::Zero();
  Energy = 0;
  
  // Se confirma que se este entre el numero maximo de iteraciones
  while (Iter<Maxit)
    {
      Iter += 1;
      //cout << "Iteracion " << Iter << endl;
      //----- Paso 2. Calcular la matriz de Fock -----//
      //--- Formamos la parte bielectrnica con P ---//
      G = Matrix2d::Zero(); //Definimos la matriz bielectronica G
      for(int i=0;i<2;i++){
	for(int j=0;j<2;j++){
	  for(int k=0;k<2;k++){
	    for(int l=0;l<2;l++){
	      G(i,j) = G(i,j)+P(l,k)*(TT[i][j][k][l]-0.5*TT[i][l][k][j]);
	    }
	  }
	}
      }
      
      
      //Sumamos HCore y G para obtener la matriz de Fock
      F = H+G;
      //cout << "P\n" << P << endl;
      //cout << "G\n" << G << endl;
      //cout << "H\n" << H << endl;
      //cout << "F\n" << F << endl;
      
      //Se calcula la energia electronica
      Energy = (P.array()*(H+F).array()*0.5).sum();
      
      //cout << "Electronic energy = " << Energy << endl;
      //----- Paso 3. Calculamos F (recuerde S^-1/2 es X y S^-1/2 es X.transpose()) -----//
      G = F*Xcanon;
      Fp = Xcanon.transpose()*G;
      
      //----- Paso 4. Diagonalizacion para encontrar el autovalor -----//
      DiagFp(Fp);
      
      //----- Paso 5. Calclamos los coeficientes de los orbitales moleculares -----//
      //Transforme eigenvectores para obtener matriz C
      C = Xcanon*Cp;
      
      //----- Paso 6. Calculamos la nueva matriz de densidad P -----//
      Oldp = P; //Guardamos la anterior matriz P
      P = Matrix2d::Zero();

      for(int i=0;i<2;i++){
	for(int j=0;j<2;j++){
	  for(int k=0;k<1;k++){
	    P(i,j) += 2.0*C(i,k)*C(j,k); //k va hasta 1 porque vamos hasta la mitad de #e
	  }
	}
      }
      
      //----- Paso 7. Verificamos el criterio de convergencia -----//
      //cout << P << endl;
      Delta = 0;
      //Se hace una suma de cuadrados y se divide por el numero de elementos en la matriz P
      Matrix2d MDelta = P-Oldp;
      Delta = sqrt(MDelta.array().pow(2).matrix().sum()/4.0);
      //cout << "Delta " << Delta << endl;
      //cout << endl;
      if(Delta<Crit)
	{
	  Energytot = Energy + Za*Zb/R; //Sumamos Eelec y Enuc
	  //cout << endl;
	  //cout << "Calculation converged with electronic energy: " << Energy << " Hartrees" << endl;
	  //cout << "Calculation converged with total energy: " << Energytot << " Hartrees" << endl;
	  //cout << endl;
	  break;
	}
    }
  
}

void HFCALC(int N,double R, double Zeta1, double Zeta2, double Za, double Zb, Matrix2d G)
{
  //Calcula las integrales
  Intgrl(N,R,Zeta1,Zeta2,Za,Zb);
  //Crea las matrices de integrales
  Colect(N,R,Zeta1,Zeta2,Za,Zb);
  //Calculo SCF
  SCF(N,R,Zeta1,Zeta2,Za,Zb,G);
}
double EMP2 = 0;
double value1 = 0, value2 = 0;
void MP2(Matrix2d C,double TT[2][2][2][2],Matrix2d E)
{
  //Me da la energia de correlacion de MP2
  double TTMO[2][2][2][2];Init2d(TTMO);
    for(int i=0;i<2;i++){
      for(int j=0;j<2;j++){
	for(int k=0;k<2;k++){
	  for(int l=0;l<2;l++){
	    for(int m=0;m<2;m++){
	      for(int n=0;n<2;n++){
		for(int o=0;o<2;o++){
		  for(int p=0;p<2;p++){
		    TTMO[i][j][k][l] += C(m,i)*C(n,j)*C(o,k)*C(p,l)*TT[m][n][o][p];
		    // cout << TTMO[i][j][k][l] << "\n";
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  int Nelec = 2;

  double integs[4][4][4][4];Init4d(integs);
  
  for(int p=1;p<5;p++){
    for(int q=1;q<5;q++){
      for(int r=1;r<5;r++){
	for(int s=1;s<5;s++){
	  value1 = TTMO[(p+1)/2-1][(r+1)/2-1][(q+1)/2-1][(s+1)/2-1] * (p%2 == r%2) * (q%2 == s%2);
	  value2 = TTMO[(p+1)/2-1][(s+1)/2-1][(q+1)/2-1][(r+1)/2-1] * (p%2 == s%2) * (q%2 == r%2);
	  integs[p-1][q-1][r-1][s-1] = value1 - value2;
	  //cout << value1 << " " << value2 << "\n";
	}
      }
    }
  }

  double fs[4];
  Vector2d dE;
  
  for(int i=0;i<4;i++)
    {
      dE = E.diagonal();
      fs[i] = dE(i/2);
    }
  //cout << fs[0] << " "  << fs[1] << " "  << fs[2] << " "  << fs[3] << " " << endl;
  EMP2 = 0.0;
  for(int i=0;i<Nelec;i++){
    for(int j=0;j<Nelec;j++){
      for(int a=Nelec;a<4;a++){
	for(int b=Nelec;b<4;b++){
	  EMP2 += 0.25*integs[i][j][a][b]*integs[i][j][a][b]/(fs[i]+fs[j]-fs[a]-fs[b]);
	}
      }
    }
  }
  
}
int main()
{
  cout.precision(8);
  //cout << scientific;
  ofstream PrintH2("data_H2.dat"),
    PrintHeH("data_HeH.dat");

  MatrixXd dataH2(81,3), dataHeH(81,3);
  
  double R = 0;
  int N = 3;
  double Zeta1 = 1.24;
  double Zeta2 = 1.24;
  double Za = 1.0;
  double Zb = 1.0;
  
  for(int n=0;n<81;n++)
    {
      R = 0.3+n*0.1;

      HFCALC(N,R,Zeta1,Zeta2,Za,Zb,G);
      HFx_axis(n) = R;
      HFener(n) = Energytot; //Matriz de energias HF
      MP2(C,TT,E);
      MP2x_axis(n) = R;
      MP2corr(n) = EMP2;
    }
  
  MP2ener = HFener + MP2corr; //Matriz de energia final MP2
  
  cout << "VALORES PARA EL H2:" << endl;
  cout << "Puntos R:\n" << "[ " << HFx_axis.transpose() << " ]" << endl;
  cout << "Energia HF:\n" << "[ " << HFener.transpose() << " ]" << endl;
  cout << "Energía correccion MP2:\n" << "[ " << MP2corr.transpose() << " ]" << endl;

  dataH2.col(0) = HFx_axis;
  dataH2.col(1) = HFener;
  dataH2.col(2) = MP2ener;
  PrintH2 << dataH2 << endl;
  
  Zeta1 = 2.0925;
  Zeta2 = 1.24;
  Za = 2.0;
  Zb = 1.0;
  
 for(int n=0;n<81;n++)
    {
      R = 0.3+n*0.1;

      HFCALC(N,R,Zeta1,Zeta2,Za,Zb,G);
      HFx_axis(n) = R;
      HFener(n) = Energytot; //Matriz de energias HF
      MP2(C,TT,E);
      MP2x_axis(n) = R;
      MP2corr(n) = EMP2;
    }
  
  MP2ener = HFener + MP2corr; //Matriz de energia final MP2
  
  cout << "VALORES PARA EL HeH+:" << endl;
  cout << "Puntos R:\n" << "[ " << HFx_axis.transpose() << " ]" << endl;
  cout << "Energia HF:\n" << "[ " << HFener.transpose() << " ]" << endl;
  cout << "Energía correccion MP2:\n" << "[ " << MP2corr.transpose() << " ]" << endl;

  dataHeH.col(0) = HFx_axis;
  dataHeH.col(1) = HFener;
  dataHeH.col(2) = MP2ener;
  PrintHeH << dataHeH << endl;

  PrintH2.close();
  PrintHeH.close();
  return 0;
}
  
