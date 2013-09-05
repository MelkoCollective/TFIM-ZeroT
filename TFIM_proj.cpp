// 
// TFIM_proj.cpp: the main body of the code
// A program to calculate the Renyi entropy in the transverse field Ising model
// Roger Melko, October 2011, University of Vermont
//
#include "head_proj.h"
#include "simparam.h"
#include "measure.h"
#include "basis.h"


int main(){


    PARAMS param; //read parameter file
    //param.printBst();
    //return 0;

    MTRand mrand(param.SEED_); //random number for metropolis

    Basis Proj(mrand);
    //Proj.printBasis();

    int Loopsize2; //loopsize squared for M^2 measurement

	if (param.EQL_ == 0)
		Proj.fileread(0);
	else{
		for (int i=0; i<param.EQL_; i++){
			Proj.DiagonalUpdate(mrand);
			//Proj.printBasis();
			Proj.LinkedList();
			//Proj.printLinkedList();
			Proj.ClusterUpdate(mrand,Loopsize2);
		}
	}

    //vector<int> inA(param.numSpin/alpha-1,0); //size of physical spins

    Measure observ;

    for (int j=0; j< param.nBin_; j++){
        observ.zero();
        for (int i=0; i<param.MCS_; i++){
            Proj.DiagonalUpdate(mrand);
            Proj.LinkedList();
            //Proj.printBasis();
            //Proj.printLinkedList();
            Proj.ClusterUpdate(mrand,Loopsize2);

            //--regular observables
            observ.measure_E(Proj);
            //observ.measure_M(Proj, Loopsize2);
            //need to do measure_M_mod before Renyi 
            observ.measure_M_mod(Proj.LeftinClust,Proj.RightinClust); //*DON'T DELETE*
            //---measure swap
            observ.Calc_Renyi(Proj.LeftinClust,Proj.RightinClust);
            //inA.assign(param.numSpin/alpha-1,0);
            //for(int k=0; k<inA.size(); k++){
            //    inA[k] = 1;
            //    //observ.Renyi_direct(k,Proj.SWAP(inA),Proj.ClustNumber);
            //    observ.Renyi_LRclust(inA,Proj.LeftinClust,Proj.RightinClust); 
            //}
            //-----------------
        }
        observ.output();
		Proj.filewrite(0); //write configuration file
    }

    return 0;

}

