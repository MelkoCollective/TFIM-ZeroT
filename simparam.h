#ifndef PARAMSIM_H
#define PARAMSIM_H

#include <fstream>
#include "head_proj.h"

//Class to read in the simulation parameters from a file
class PARAMS
{
    public:
        int numSpin; 
        int numLattB; //number of lattice bonds is 2N
		int EQL_;     //the number of Monte Carlo steps
		int MCS_;     //the number of Monte Carlo steps
        int nX_;      //linear size of lattice
        int nY_;      //linear size of lattice
        int nBin_;    //number of production bins
        long SEED_;

        long m_;     //The order of the expansion
        double h_x;  //Transverse field 

        PARAMS();

        //lattice spin coordination numbers
        vector <index2> Bst;

        void printBst();

		//stuff for the swap
		int ratioON;  //1 for ratio, 0 for bare swap, from regionX.dat file
		int nSwap;
        int numRealSpin; //the number of spins in the physical system
		vector<vector<int> > inAreg; //definition of the A regions to measure

   
}; //PARAMS

PARAMS::PARAMS(){
    //initializes commonly used parameters from a file
    ifstream pfin;
    pfin.open("param.dat");

    pfin >> nX_;
    pfin >> nY_;
    pfin >> EQL_;
    pfin >> MCS_;
    pfin >> nBin_;
    pfin >> SEED_;
    pfin >> m_;
    pfin >> h_x;
    pfin.close();

    if (nY_ == 0) {    // ---------1D chain
        //derived constants
        numSpin = alpha*nX_;    //replicate it: alpha is a global #define
        numLattB = numSpin; //Periodic BC

        //Initialize lattice bond array
        int a,b;
        index2 temp;
        for (int rep=0; rep<alpha; rep++)
            for (int i=0; i<numLattB/alpha; i++){  //EACH LAYER
                a = rep*numLattB/alpha+i;
                b = a+1;
                if (b == (rep+1)*numSpin/alpha) b = rep*numLattB/alpha; //PBC
                temp.set(a,b);
                Bst.push_back(temp);
            }//i

    }//1D chain

    else{//2D

        if (alpha > 2) cout<<"S_alpha greater than S_2 NOT TESTED FOR 2D \n";
        if (nX_ != nY_) cout<<"Rectangular lattices not tested \n";
        numSpin = alpha*nX_*nX_;
        numLattB = 2*numSpin; //Periodic BC for 2D lattice

        //Initialize lattice bond array
        int a,b,d;
        index2 temp;
        for (int rep=0; rep<alpha; rep++)
            for (int i=0; i<numSpin/alpha; i++){  
                //horizontal bond
                a = rep*numSpin/alpha + i;
                b = a+1;
                if ( b%nX_ == 0) b -= nX_;
                temp.set(a,b);
                Bst.push_back(temp);
                //vertical bond
                a = rep*numSpin/alpha + i;
                d = a+nX_;
                if (d>= rep*numSpin/alpha + nX_*nX_) d -= nX_*nX_;
                temp.set(a,d);
                Bst.push_back(temp);
            }//i

    }//end 2D

	//--------read in regions A and X
	vector<int> Atemp;  //vector to be pushed back
	Atemp.assign(numSpin/alpha,0);

	ifstream fin;
	fin.open("regionA.dat");

	if (fin.fail() ) { //check for errors
		cout<<"Could not open a regionA.dat file"<<endl;
	}

	fin>>nSwap;
	if (nSwap < 1) cout<<"regionA.dat error 1 \n";

	int temp;

	for (int j=0; j<nSwap; j++){
		for (int i=0; i<numSpin/alpha; i++){
			fin>>temp;
			if (temp != 0 && temp != 1)  cout<<"regionA.dat error 2 \n";
			Atemp.at(i) = temp; 
		}
		fin>>temp;
		if (temp != -99) cout<<"regionA.dat error 3 \n";
		inAreg.push_back(Atemp); //the 1 spin region
	}//j

	fin.close();

	//--if there is a regionX.dat file, push this back as the last
	//--vector element in inAreg<>
    fin.open("regionX.dat");
    if (fin.fail() ) { //check for errors
        cout<<"Could not open a regionX.dat file"<<endl;
        ratioON = 0;
    }
    else{
        fin>>temp;
        if ( temp==0 ){  //check for errors
            //cout<<"not using Ratio! : renyi"<<endl;
            ratioON = 0;
        }
        else{
            ratioON = 1;
            for (int i=0; i<numSpin/alpha; i++){
                fin>>temp;
                if (temp != 0 && temp != 1)  cout<<"regionX.dat error 3  \n";
                Atemp.at(i) = temp; 
            }
            fin>>temp;
            if (temp != -99) cout<<"regionX.dat error 4  BASIS\n";
            inAreg.push_back(Atemp); //the 1 spin region
        }

        fin.close();
    }//regionX.dat found

	//------ done reading in regions A and X

    //for (int i=0; i<inAreg.size(); i++){
    //    for (int j=0; j<inAreg[0].size(); j++)
    //        cout<<inAreg[i][j]<<" ";
    //    cout<<endl;
    //}

    numRealSpin = numSpin/alpha;

}//constructor

void PARAMS::printBst(){
    cout<<"Bst "<<endl;
    for (int i=0; i<Bst.size(); i++)
        Bst[i].print();
}

#endif
