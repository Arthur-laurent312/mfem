
#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

int main(int argc, char *argv[])
{
	ifstream ex2_boucle("sol_boucle.txt");
	ifstream ex2("sol.txt");
   if(ex2_boucle)
   {
	if(ex2)
{
	string ligne;
	//Taille fichier 1
	int size = 0;
	while(getline(ex2_boucle, ligne))
	{
		size++;
	}
	Vector sol_boucle(size);
	ex2_boucle.clear();
	ex2_boucle.seekg( 0,ios::beg);

	 size = 0;
      while(getline(ex2_boucle, ligne))
      {
		sol_boucle(size) = stod(ligne);
		size++;
      }
	ex2_boucle.clear();
	ex2_boucle.seekg( 0,ios::beg);
	//Taille fichier 2
	 int size1 = 0;
      while(getline(ex2_boucle, ligne))
      {
		size1++;
      }

	Vector sol(size1);
	ex2.clear();
	ex2.seekg( 0,ios::beg);

	size1=0;
      while(getline(ex2, ligne))
      {
		sol(size1) = stod(ligne);
		size1++;

      }

	if(size1==size){
		sol -= sol_boucle;
		cout<<sol.Norml2()<<endl;
	}
	else{
		cout<<"Problème de taille de fichier"<<endl;
	}

   }
   else
   {
      cout << "ERREUR: Impossible d'ouvrir le fichier en lecture." << endl;
   }

   return 0;
}
}

