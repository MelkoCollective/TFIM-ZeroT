#ifndef BASIS_H
#define BASIS_H

#include "head_proj.h"
#include "simparam.h"


class Basis: public PARAMS
{
    private:
       vector<int> LinkList;
       vector<int> LegType; 
       vector<index3> Associates;
	   vector<int> inCluster;

	   vector<int> CalcRightCluster();

	public:
	   vector <index2> OperatorList; //The operator list of 2m
	   // (-2,i): an off-diagonal site operator h(sigma^+_i + sigma^-_i)
	   // (-1,i): a diagonal site operator h
	   // (i,j):  a diagonal bond operator J(sigma^z_i sigma^z_j + 1)

	   vector<int> S_left; //the left and right trial spin state
	   vector<int> S_right;

	   //center clusters built in Linked List
	   vector<int> LeftinClust;
	   vector<int> RightinClust;

	   int ClustNumber; //the number of clusters in the center

	   Basis(MTRand &); //constructor
	   void DiagonalUpdate(MTRand &);
	   void LinkedList();
	   void ClusterUpdate(MTRand &, int&);
	   int calc_LoopSize2();
	   int SWAP(const vector<int>& );
	   void printBasis();
	   void printLinkedList();
	   void filewrite(const int & num);
	   void fileread(const int & num);


};


Basis::Basis(MTRand& ran){//constructor

    for (int i=0; i< numSpin; i++){
        S_left.push_back(1);   //start with all spin parallel
        S_right.push_back(1);
    }

    int bond;
    index2 temp;
    for (int i=0; i<2*m_; i++){

        if (ran.randInt(1) == 0) { //flip a coin
            temp.set(-1,ran.randInt(numSpin-1)); //diagonal site operator
        }
        else {
            bond = ran.randInt(numLattB-1);      //diagonal bond operator
            temp = Bst[bond]; //checked overloaded =
            //cout<<"temp "<<bond<<" ";
            //Bst[bond].print();
            //cout<<endl;
        }

        OperatorList.push_back(temp);

    }//i
     
    //S_right.at(S_right.size()-1) = S_right.at(S_right.size()-1)^1; 
    //temp.set(-2,numSpin-1);
    //OperatorList.at(m_) =  temp;

}//----------------constructor



//----------------DiagonalUpdate function
void Basis::DiagonalUpdate(MTRand& ran){

    vector<int> S_prop; //this is the temporary propagated spin state
    S_prop = S_left; //assign to the left spin state

    //Probability for a single-site diagonal operator
    //double hProb = h_x/(2.0+h_x);
    double hProb = 1.0*numSpin*h_x/(numLattB*2.0+numSpin*h_x);

    int bond;
    int flag;
    for(int i=0; i<OperatorList.size(); i++){

        if (ratioON == 1 && i == OperatorList.size()/2){ //SWAP for ratio trick
            int tempj;
            int Xindex = inAreg.size()-1;
            for (int j=0; j<inAreg[Xindex].size(); j++)
                if (inAreg[Xindex][j] == 1){
                    tempj = S_prop[j];
                    for (int rep=0; rep<alpha-1; rep++) //Permutation
                        S_prop[rep*numRealSpin+j] = S_prop[rep*numRealSpin+j+numRealSpin];
                    S_prop[(alpha-1)*numRealSpin+j] = tempj;
                    //S_prop[j] = S_prop[j+numSpin/alpha]; //old SWAP
                    //S_prop[j+numSpin/alpha] = tempj;
                }
        }//--- ratio trick

        if (OperatorList[i].A == -2) //this is a off-diagonal site operator
            S_prop[OperatorList[i].B] = S_prop[OperatorList[i].B]^1 ; //binary spin flip

        else { //sample a new diagonal operator

            flag = 0;
            do{ //repeat until a valid selection is made
                if (hProb > ran.rand() ){ //probability to choose a single-site h operator
                    OperatorList[i].set(-1,ran.randInt(numSpin-1)); //diagonal site operator
                    flag = 1; //successful!
                }
                else{
                    //cout<<"B "<<i<<" "<<endl;
                    bond = ran.randInt(numLattB-1);  //new bond for diagonal bond operator
                    if (S_prop[Bst[bond].A] == S_prop[Bst[bond].B]) {
                        //spins are the same on the new bond
                        OperatorList[i] = Bst[bond]; //check overloaded =
                        flag = 1; //successful!
                    }
                } //else the diagonal operator stays unchanged!  repeat
            }while(flag == 0);
        }//insert diagonal

    }//i  the 2*m propagation

    //DEBUG: check if the state was propagated correctly
    for (int i=0; i<S_prop.size(); i++)
        if (S_prop[i] != S_right[i]) cout<<"Basis state prop error: DIAG UPDATE \n";


}//----------------DiagonalUpdate



//----------------LinkedList function
void Basis::LinkedList(){

    //Clear linked list from last iteration
    LinkList.clear(); 
    LegType.clear();
    Associates.clear();
    
    index3 empty(-1,-1,-1);
    index3 temp3;

    vector<int> First;
    for (int i=0; i<numSpin; i++){ //the first vertex leg for each spin
        First.push_back(i);
        //below, build the first N vertices from the left-basis
        LinkList.push_back(-99); //unknown what these link to!
        LegType.push_back(S_left[i]); //0 or 1
        Associates.push_back(empty); //these have no associates
    }

    vector<int> S_prop; //this is the temporary propagated spin state
    S_prop = S_left; //assign to the left spin state

    int Ccount, Ctemp; //cluster counter
	vector<int> LRinClust;
	Ccount = 0;
	LRinClust.assign(numSpin,0);

    int count = numSpin;  
    int site, site1, site2;
    //The linked list is now size N.  Add the 2m operators each of 4 or 2 legs
    for(int i=0; i<OperatorList.size(); i++){

		if (ratioON == 1 && i == OperatorList.size()/2){ //SWAP for ratio trick
				int tempj;
				int Xindex = inAreg.size()-1;
				for (int j=0; j<inAreg[Xindex].size(); j++)
					if (inAreg[Xindex][j] == 1){
                        tempj = First[j];
                        for (int rep=0; rep<alpha-1; rep++) //Permutation
                            First[rep*numRealSpin+j] = First[rep*numRealSpin+j+numRealSpin];
                        First[(alpha-1)*numRealSpin+j] = tempj;
						//First[j] = First[j+numSpin/alpha]; //old SWAP
						//First[j+numSpin/alpha] = tempj;
                        tempj = S_prop[j];
                        for (int rep=0; rep<alpha-1; rep++) //Permutation
                            S_prop[rep*numRealSpin+j] = S_prop[rep*numRealSpin+j+numRealSpin];
                        S_prop[(alpha-1)*numRealSpin+j] = tempj;
                        //S_prop[j] = S_prop[j+numSpin/alpha]; //old SWAP
                        //S_prop[j+numSpin/alpha] = tempj;
					}
		}//ratio trick
       
        //assign non-trivial associates
        if (OperatorList[i].A != -2 && OperatorList[i].A != -1){
            temp3.set(count+1,count+2,count+3); Associates.push_back(temp3);
            temp3.set(count,count+2,count+3); Associates.push_back(temp3);
            temp3.set(count,count+1,count+3); Associates.push_back(temp3);
            temp3.set(count,count+1,count+2); Associates.push_back(temp3);
            count += 4;
        }//done assigning associates
        else{
            Associates.push_back(empty);
            Associates.push_back(empty);
            count += 2;
        }

        if (OperatorList[i].A == -2){ //1-site off-diagonal operator is encountered
            site = OperatorList[i].B;
            //"lower" or leftmost leg
            LinkList.push_back(First[site]); //site index
            LegType.push_back(S_prop[site]); //the spin of the leg
            S_prop[site] = S_prop[site]^1;   //this is off-d: flip it
            LinkList[First[site]] = LinkList.size()-1; //this leg links backwards...
            First[site] = LinkList.size(); //update
            //"upper" or rightmost leg
            LinkList.push_back(-99); //null site index
            LegType.push_back(S_prop[site]); //the spin of the leg (flipped)
            LRinClust[site] = 0; //cluster counter for LHS
        }
        else if (OperatorList[i].A == -1){ //1-site diagonal operator is encountered
            site = OperatorList[i].B;
            //"lower" or leftmost leg
            LinkList.push_back(First[site]); //site index
            LegType.push_back(S_prop[site]); //the spin of the leg
            LinkList[First[site]] = LinkList.size()-1; //this leg links backwards...
            First[site] = LinkList.size(); //update
            //"upper" or rightmost leg
            LinkList.push_back(-99); //null site index
            LegType.push_back(S_prop[site]); //the spin of the leg 
            LRinClust[site] = 0; //cluster counter for LHS
        }
        else {//2-site diagonal operator is encountered (4 legs)
            //lower left
            site1 = OperatorList[i].A;
            LinkList.push_back(First[site1]); //site index
            LegType.push_back(S_prop[site1]); //the spin of the leg
            LinkList[First[site1]] = LinkList.size()-1; //this leg links backwards...
            First[site1] = LinkList.size()+1;
            //lower right
            site2 = OperatorList[i].B;
            LinkList.push_back(First[site2]); //site index
            LegType.push_back(S_prop[site2]); //the spin of the leg
            LinkList[First[site2]] = LinkList.size()-1; //this leg links backwards...
            First[site2] = LinkList.size()+1;
            //upper left
            LinkList.push_back(-99); //null site index
            LegType.push_back(S_prop[site1]); //the spin of the leg 
            //upper right
            LinkList.push_back(-99); //null site index
            LegType.push_back(S_prop[site2]); //the spin of the leg 
            //------Build the clusters here
            //cout<<"sites "<<site1<<" "<<site2<<endl;
            //cout<<"clust "<<LRinClust[site1]<<" "<<LRinClust[site2]<<endl;
            if (LRinClust[site1] == 0 && LRinClust[site2] == 0)
            {
                Ccount++;
                LRinClust[site1] = Ccount; LRinClust[site2] = Ccount;
            }
            else if (LRinClust[site1] != 0 && LRinClust[site2] == 0)
                LRinClust[site2] = LRinClust[site1];
			else if (LRinClust[site2] != 0 && LRinClust[site1] == 0)
				LRinClust[site1] = LRinClust[site2];
			else if (LRinClust[site2] != 0 && LRinClust[site1] != 0)
			{
				if (LRinClust[site2] != LRinClust[site1] ){
					Ctemp = LRinClust[site2];
					for (int ii=0; ii<LRinClust.size(); ii++)
						if (LRinClust[ii] == Ctemp )
							LRinClust[ii] = LRinClust[site1];
				}
			}
            else cout<<"Mid cluster error \n";
            //cout<<"Aclust "<<LRinClust[site1]<<" "<<LRinClust[site2]<<endl;
            //-------------------------
        }//bond operator encountered

        if (i == OperatorList.size()/2-1){
            LeftinClust = LRinClust; //Left-hand side
            for (int jj=0; jj<LeftinClust.size(); jj++){
                if (LeftinClust[jj] == 0){
                    Ccount++;
                    LeftinClust[jj] = Ccount;
                }
            }//jj
        }

    }//i

    //RightinClust = LRinClust;  //Right-hand side
    //for (int jj=0; jj<RightinClust.size(); jj++){
    //    if (RightinClust[jj] == 0 || RightinClust[jj] == -1){
    //        Ccount++;
    //        RightinClust[jj] = Ccount;
    //    }
    //}//jj

	RightinClust = CalcRightCluster();

    //now add the legs of the final ("top"or right-hand) spin state
    for (int i=0; i<numSpin; i++){ 
        LinkList.push_back(First[i]);
        LinkList[First[i]] = LinkList.size()-1;
        LegType.push_back(S_prop[i]); //0 or 1
        Associates.push_back(empty);
    }

//    cout<<"Ass size :"<<Associates.size()<<endl;
//    cout<<"LL size :"<<LinkList.size()<<endl;
//    cout<<"LT size :"<<LegType.size()<<endl;

    //DEBUG: check if the state was propagated correctly
    for (int i=0; i<S_prop.size(); i++)
        if (S_prop[i] != S_right[i]) cout<<"Basis state prop error: LINKED LIST\n";

}//----------------LinkedList

vector<int> Basis::CalcRightCluster(){

    vector<int> Last;
    vector<int> LL;

    for (int i=0; i<numSpin; i++){ //the first vertex leg for each spin
        Last.push_back(i);
        LL.push_back(-99); //unknown what these link to!
	}

    int Ccount = 0;
	int Ctemp;
	vector<int> LRinClust(numSpin,0);

    int site, site1, site2;
    for(int i=OperatorList.size()-1; i>=OperatorList.size()/2; i--){

        if (OperatorList[i].A == -2){ //1-site off-diagonal operator is encountered
            site = OperatorList[i].B;
            //"lower" or leftmost leg
            LL.push_back(Last[site]); //site index
            LL[Last[site]] = LL.size()-1; //this leg links backwards...
            Last[site] = LL.size(); //update
            //"upper" or rightmost leg
			LL.push_back(-99); //null site index
			LRinClust[site] = 0; //cluster counter for LHS
        }
        else if (OperatorList[i].A == -1){ //1-site diagonal operator is encountered
            site = OperatorList[i].B;
            //"lower" or leftmost leg
            LL.push_back(Last[site]); //site index
            LL[Last[site]] = LL.size()-1; //this leg links backwards...
            Last[site] = LL.size(); //update
            //"upper" or rightmost leg
            LL.push_back(-99); //null site index
		    LRinClust[site] = 0; //cluster counter for LHS

        }
        else {//2-site diagonal operator is encountered (4 legs)
            //lower left
            site1 = OperatorList[i].A;
            LL.push_back(Last[site1]); //site index
            LL[Last[site1]] = LL.size()-1; //this leg links backwards...
            Last[site1] = LL.size()+1;
            //lower right
            site2 = OperatorList[i].B;
            LL.push_back(Last[site2]); //site index
            LL[Last[site2]] = LL.size()-1; //this leg links backwards...
            Last[site2] = LL.size()+1;
            //upper left
            LL.push_back(-99); //null site index
            //upper right
			LL.push_back(-99); //null site index
			if (LRinClust[site1] == 0 && LRinClust[site2] == 0)
			{
				Ccount++;
				LRinClust[site1] = Ccount; LRinClust[site2] = Ccount;
			}
			else if (LRinClust[site1] != 0 && LRinClust[site2] == 0)
				LRinClust[site2] = LRinClust[site1];
			else if (LRinClust[site2] != 0 && LRinClust[site1] == 0)
				LRinClust[site1] = LRinClust[site2];
			else if (LRinClust[site2] != 0 && LRinClust[site1] != 0)
			{
				if (LRinClust[site2] != LRinClust[site1] ){
					Ctemp = LRinClust[site2];
					for (int ii=0; ii<LRinClust.size(); ii++)
						if (LRinClust[ii] == Ctemp )
							LRinClust[ii] = LRinClust[site1];
				}
			}
			else cout<<"Mid cluster error \n";
        }//bond operator encountered
    }//i

	for (int jj=0; jj<LRinClust.size(); jj++){
		if (LRinClust[jj] == 0){
			Ccount++;
			LRinClust[jj] = Ccount;
		}
	}//jj

    //A permute for the ratio trick...
    if (ratioON == 1){
        int Xindex = inAreg.size()-1;
        int temp;
        for (int i=0; i<inAreg[Xindex].size(); i++){
            if (inAreg[Xindex][i] != 0){ 
                temp = LRinClust[i];
                for (int rep=0; rep<alpha-1; rep++)
                    LRinClust[rep*numRealSpin+i] = LRinClust[rep*numRealSpin+i+numRealSpin];
                LRinClust[(alpha-1)*numRealSpin+i] = temp;
            }
        }
    }//if

	return LRinClust;


}//CalcRightCluster

//----------------ClusterUpdate
void Basis::ClusterUpdate(MTRand& ran, int& L2){

    inCluster.clear();  //redundant with assign
    inCluster.assign(LinkList.size(),0);//nothing in clusters yet

    stack<int> cluster;

    int leg, assoc;
    bool flip;

    int ccount = 0;
    for (int i=0; i<LinkList.size(); i++){ //loop to find all clusters

        //add a new leg
        if (inCluster[i] == 0 && Associates[i].A == -1){ //spins and site ops only
            ccount ++; //cluster counter

            cluster.push(i);
            inCluster[cluster.top()] = ccount;

            if (ran.rand() < 0.5) flip = true; else flip = false; //flip a coin for SW

            if (flip == true) LegType[cluster.top()] = LegType[cluster.top()]^1;

            while(!cluster.empty()){ //build the cluster associated with this leg

                //first follow the link
                leg = LinkList[cluster.top()];
                cluster.pop();

                if (inCluster[leg] == 0){
                    inCluster[leg] = ccount; //add the linked leg
                    if (flip == true) LegType[leg] = LegType[leg]^1;
                    //now check all associates
                    assoc = Associates[leg].A;
                    if (assoc != -1) { 
                        cluster.push(assoc); inCluster[assoc] = ccount; 
                        if (flip == true) LegType[assoc] = LegType[assoc]^1;
                        assoc = Associates[leg].B;
                        cluster.push(assoc); inCluster[assoc] = ccount; 
                        if (flip == true) LegType[assoc] = LegType[assoc]^1;
                        assoc = Associates[leg].C;
                        cluster.push(assoc); inCluster[assoc] = ccount; 
                        if (flip == true) LegType[assoc] = LegType[assoc]^1;
                    }
                }

            }//while

        }//if building a new cluster

    }//i

    //cout<<"inCluster: ";
    //for (int i=0; i<inCluster.size(); i++)
    //    cout<<inCluster[i]<<" ";
    //cout<<endl;

    //map back basis states and operator list
    for(int i=0; i<numSpin; i++){
        S_left[i] = LegType[i];
    }//i
    int count = numSpin;  
    for(int i=0; i<OperatorList.size(); i++){

        //assign non-trivial associates
        if (OperatorList[i].A != -2 && OperatorList[i].A != -1){
            count += 4;
        }//done assigning associates
        else{
            if (LegType[count] == LegType[count+1])
                OperatorList[i].A = -1;
            else
                OperatorList[i].A = -2;
            count += 2;
        }
    }//i
    for(int i=0; i<numSpin; i++){
        S_right[i] = LegType[LegType.size()-numSpin + i];
    }//i

    //L2 = calc_LoopSize2(); //used for old magnetization estimator

    //LinkList.clear(); //clear up the linked list
    //LegType.clear();
    //Associates.clear();
    //inCluster.clear();  //redundant with assign

}//----------------ClusterUpdate

int Basis::calc_LoopSize2(){

    vector<int> Last;
    vector<int> LL;

    for (int i=0; i<numSpin; i++){ //the first vertex leg for each spin
        Last.push_back(i);
        LL.push_back(-99); //unknown what these link to!
    }

    int site, site1, site2;
    for(int i=0; i<OperatorList.size()/2; i++){

        if (OperatorList[i].A == -2){ //1-site off-diagonal operator is encountered
            site = OperatorList[i].B;
            //"lower" or leftmost leg
            LL.push_back(Last[site]); //site index
            LL[Last[site]] = LL.size()-1; //this leg links backwards...
            Last[site] = LL.size(); //update
            //"upper" or rightmost leg
            LL.push_back(-99); //null site index
        }
        else if (OperatorList[i].A == -1){ //1-site diagonal operator is encountered
            site = OperatorList[i].B;
            //"lower" or leftmost leg
            LL.push_back(Last[site]); //site index
            LL[Last[site]] = LL.size()-1; //this leg links backwards...
            Last[site] = LL.size(); //update
            //"upper" or rightmost leg
            LL.push_back(-99); //null site index
        }
        else {//2-site diagonal operator is encountered (4 legs)
            //lower left
            site1 = OperatorList[i].A;
            LL.push_back(Last[site1]); //site index
            LL[Last[site1]] = LL.size()-1; //this leg links backwards...
            Last[site1] = LL.size()+1;
            //lower right
            site2 = OperatorList[i].B;
            LL.push_back(Last[site2]); //site index
            LL[Last[site2]] = LL.size()-1; //this leg links backwards...
            Last[site2] = LL.size()+1;
            //upper left
            LL.push_back(-99); //null site index
            //upper right
            LL.push_back(-99); //null site index
        }
    }//i

    vector<int> ClusterSize(LL.size(),0); //is this the maximum size?
    vector<int> numCluster(LL.size(),0); 

    for (int i=0; i<numSpin; i++){
        ClusterSize[inCluster[Last[i]]]++;
        numCluster[inCluster[Last[i]]]=1;
    }

    //cout<<"ClusterSize: ";
    //for (int i=0; i<ClusterSize.size(); i++)
    //    if (ClusterSize[i] != 0) cout<<i<<" "<<ClusterSize[i]<<" ";
    //cout<<endl;

    int sizesquared=0;
    int Ccount=0;
    for (int i=0; i<ClusterSize.size(); i++){
        sizesquared += ClusterSize[i]*ClusterSize[i];
        Ccount += numCluster[i];
    }

    //cout<<sizesquared<<endl;
    //cout<<"clusters in center "<<Ccount<<" ";

    ClustNumber = Ccount; //number of clusters crossing the center 
    return sizesquared;

}//calc_LoopSize2


//----------------print LinkedList
void Basis::printLinkedList(){

    //for (int i=0; i<LinkList.size(); i++){ 
    //    cout<<i<<" ";
    //    cout<<LinkList[i]<<" ";
    //    cout<<LegType[i]<<"\n";
    //}

    cout<<"LH cluster : ";
    for (int i=0; i<LeftinClust.size(); i++)
        cout<<LeftinClust[i]<<" ";
    cout<<endl;

    cout<<"RH cluster : ";
    for (int i=0; i<RightinClust.size(); i++)
        cout<<RightinClust[i]<<" ";
    cout<<endl;

}//printLinkedList



//----------------print function
void Basis::printBasis(){
    cout<<"Basis "<<endl;
    for (int i=0; i<S_left.size(); i++){
        cout<<S_left[i]<<" ";
    }
    cout<<endl;
    for (int i=0; i<S_right.size(); i++){
        cout<<S_right[i]<<" ";
    }
    cout<<endl;
    for (int i=0; i<OperatorList.size(); i++){
        OperatorList[i].print();
    }
//    cout<<endl;
//    for (int i=0; i<Associates.size(); i++){
//        Associates[i].print();
//    }
//    cout<<endl;

}//print


//----------------SWAP
int Basis::SWAP(const vector<int>& inA){

    vector<int> Last;
    vector<int> swapLL;
    for (int i=0; i<numSpin; i++){ //the first vertex leg for each spin
        Last.push_back(i);
        swapLL.push_back(-99);     //unknown what these link to!
    }

    vector<int> Last_atHalf;

    int count = numSpin;  
    int site, site1, site2;
    for(int i=0; i<OperatorList.size(); i++){
         
        if (i == OperatorList.size()/2){ //check - in the center?
            int tempj;
            for (int j=0; j<inA.size(); j++)
                if (inA[j] == 1){
                    tempj = Last[j];
                    Last[j] = Last[j+numSpin/2]; //SWAP: not for Permutation
                    Last[j+numSpin/2] = tempj;
                }
            Last_atHalf = Last; //copy 
        }//i at the center

        if (OperatorList[i].A == -2){ //1-site off-diagonal operator is encountered
            site = OperatorList[i].B;
            //"lower" or leftmost leg
            swapLL.push_back(Last[site]); //site index
            swapLL[Last[site]] = swapLL.size()-1; //this leg links backwards...
            Last[site] = swapLL.size(); //update
            //"upper" or rightmost leg
            swapLL.push_back(-99); //null site index
        }
        else if (OperatorList[i].A == -1){ //1-site diagonal operator is encountered
            site = OperatorList[i].B;
            //"lower" or leftmost leg
            swapLL.push_back(Last[site]); //site index
            swapLL[Last[site]] = swapLL.size()-1; //this leg links backwards...
            Last[site] = swapLL.size(); //update
            //"upper" or rightmost leg
            swapLL.push_back(-99); //null site index
        }
        else {//2-site diagonal operator is encountered (4 legs)
            //lower left
            site1 = OperatorList[i].A;
            swapLL.push_back(Last[site1]); //site index
            swapLL[Last[site1]] = swapLL.size()-1; //this leg links backwards...
            Last[site1] = swapLL.size()+1;
            //lower right
            site2 = OperatorList[i].B;
            swapLL.push_back(Last[site2]); //site index
            swapLL[Last[site2]] = swapLL.size()-1; //this leg links backwards...
            Last[site2] = swapLL.size()+1;
            //upper left
            swapLL.push_back(-99); //null site index
            //upper right
            swapLL.push_back(-99); //null site index
        }

    }//i

    //now add the legs of the final ("top"or right-hand) spin state
    for (int i=0; i<numSpin; i++){ 
        swapLL.push_back(Last[i]);
        swapLL[Last[i]] = swapLL.size()-1;
    }
    //DONE BUILDING SWAPPED LINKED LIST
    
    vector<int> swap_inClust(swapLL.size(),0);//nothing in clusters yet
    stack<int> cluster;

    int leg, assoc;
    int ccount = 0;
    for (int i=0; i<swapLL.size(); i++){ //loop to find all clusters
        //add a new leg
        if (swap_inClust[i] == 0 && Associates[i].A == -1){ //spins and site ops only
            ccount ++; //cluster counter
            cluster.push(i);
            swap_inClust[cluster.top()] = ccount;

            while(!cluster.empty()){ //build the cluster associated with this leg
                //first follow the link
                leg = swapLL[cluster.top()];
                cluster.pop();

                if (swap_inClust[leg] == 0){
                    swap_inClust[leg] = ccount; //add the linked leg
                    //now check all associates
                    assoc = Associates[leg].A;
                    if (assoc != -1) { 
                        cluster.push(assoc); swap_inClust[assoc] = ccount; 
                        assoc = Associates[leg].B;
                        cluster.push(assoc); swap_inClust[assoc] = ccount; 
                        assoc = Associates[leg].C;
                        cluster.push(assoc); swap_inClust[assoc] = ccount; 
                    }
                }
            }//while

        }//if building a new cluster
    }//i

    //DONE BUILDING CLUSTERS

    vector<int> numCluster(swapLL.size(),0); //is this the maximum size?

    for (int i=0; i<numSpin; i++)
       numCluster[swap_inClust[Last_atHalf[i]]] = 1; //look at half

    int Ccount=0;
    for (int i=0; i<numCluster.size(); i++)
        Ccount += numCluster[i];

    return Ccount;

}//----------------SWAP


void Basis::filewrite(const int & num){

    char fname[8];

    if (num == 0) fname[1] = '0';
    else if (num%9 == 0) fname[1] = '9';
    else if (num%8 == 0) fname[1] = '8';
    else if (num%7 == 0) fname[1] = '7';
    else if (num%6 == 0) fname[1] = '6';
    else if (num%5 == 0) fname[1] = '5';
    else if (num%4 == 0) fname[1] = '4';
    else if (num%3 == 0) fname[1] = '3';
    else if (num%2 == 0) fname[1] = '2';
    else if (num%1 == 0) fname[1] = '1';

    fname[0] = '0';
    fname[2] = '.';
    fname[3] = 'b';
    fname[4] = 'a';
    fname[5] = 's';
    fname[6] = 'e';
    fname[7] = '\0';

    ofstream cfout;
    cfout.open(fname);

	for (int i=0; i<S_left.size(); i++)
		cfout<<S_left[i]<<" ";
	cfout<<endl;
	for (int i=0; i<S_right.size(); i++)
		cfout<<S_right[i]<<" ";
	cfout<<endl;

	cfout<<OperatorList.size()<<endl;

	for (int i=0; i<OperatorList.size(); i++){
		cfout<<OperatorList[i].A<<" ";
		cfout<<OperatorList[i].B<<endl;
	}

	cfout<<"-999 \n"; //check for file corruption etc

    cfout.close();

}//filewrite


void Basis::fileread(const int & num){

    char fname[8];

    if (num == 0) fname[1] = '0';
    else if (num%9 == 0) fname[1] = '9';
    else if (num%8 == 0) fname[1] = '8';
    else if (num%7 == 0) fname[1] = '7';
    else if (num%6 == 0) fname[1] = '6';
    else if (num%5 == 0) fname[1] = '5';
    else if (num%4 == 0) fname[1] = '4';
    else if (num%3 == 0) fname[1] = '3';
    else if (num%2 == 0) fname[1] = '2';
    else if (num%1 == 0) fname[1] = '1';

    fname[0] = '0';
    fname[2] = '.';
    fname[3] = 'b';
	fname[4] = 'a';
	fname[5] = 's';
	fname[6] = 'e';
	fname[7] = '\0';

	ifstream cfin;
	cfin.open(fname);

	if (cfin.fail() ) { //check for errors
		cout<<"Could not open a basis input file "<<endl;
	}

	int temp;
	for (int i=0; i<S_left.size(); i++){
		cfin>>temp;
		S_left[i] = temp;
	}
	for (int i=0; i<S_right.size(); i++){
		cfin>>temp;
		S_right[i] = temp;
	}

	cfin>>temp;  //OperatorList.size();
	if (temp != 2*m_) cout<<"Basis fileread error 1 \n";
	if (temp != OperatorList.size() ) cout<<"Basis fileread error 2 \n";

	for (int i=0; i<OperatorList.size(); i++){
		cfin>>temp;
		OperatorList[i].A = temp;
		cfin>>temp;
		OperatorList[i].B = temp;
	}

	cfin>>temp;  //OperatorList.size();
	if (temp != -999) cout<<"Basis fileread error 3 \n";

	cfin.close();


}//fileread


#endif
