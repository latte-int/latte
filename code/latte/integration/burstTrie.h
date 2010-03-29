#ifndef BURSTTRIE_H
#define BURSTTRIE_H

#define BURST_MAX 10
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <stdio.h>
#include <sstream>

NTL_CLIENT

//this defines the data structures used
template <class T, class S> struct term
{
	T* coef;
	S* exponents;
	int expCount;
        int degree;
	
	term<T, S>* nextTerm;
};

template <class T, class S> struct burstTrie;
template <class T, class S> struct container
{
	//T* myCoef; //coefficient of the monomial in the sorted variable alone
        //int myDegree; //degree of the monomial in the sorted variable alone
        term<T, S>* myTerm; //use this instead of the 2 above, set expCount to 0 and don't bother initializing exponents
        
	term<T, S>* firstTerm; //head of linked list of terms, sorted on remaining variables
	burstTrie<T, S>* myTrie; //burst threshold exceeded, convert to burstTrie
	int itemCount;
	
	container<T, S>* nextCont;
};

template <class T, class S> struct burstTrie
{
	S* range; 
	container<T, S>* firstCont; //head of sorted linked list of containers (sortMax - sortMin elements)
};

//Constructors
template <class T, class S>
container<T, S>* createContainer()
{
	container<T, S>* curContainer;
	curContainer = (container<T, S>*) malloc (sizeof(container<T, S>));
	curContainer->nextCont = NULL;
	curContainer->itemCount = 0;
	curContainer->myTrie = NULL;
        curContainer->myTerm = NULL;
	return curContainer;
}

template <class T, class S>
term<T, S>* stripTerm(term<T, S>* oldTerm)
{
    if (oldTerm->expCount == 0) { cout << "Stripping coeff" << endl; }
	term<T, S>* newTerm = (term<T, S>*) malloc (sizeof(term<T, S>));
	newTerm->coef = new T(*oldTerm->coef);
	newTerm->expCount = oldTerm->expCount - 1;
	newTerm->exponents = new S[newTerm->expCount];
	for (int i = 0; i < newTerm->expCount; i++)
	{ newTerm->exponents[i] = oldTerm->exponents[i + 1]; }
        newTerm->degree = oldTerm->degree;
	return newTerm;
}

//Destructors
template <class T, class S>
void destroyTrie(burstTrie<T, S>* myTrie)
{
	cout << "Destroying trie" << endl;
	container<T, S> *curCont, *tempCont;
	curCont = myTrie->firstCont;
	for (S i = myTrie->range[0]; i <= myTrie->range[1]; i++)
	{
		tempCont = curCont;
		curCont = curCont->nextCont;
		destroyContainer(tempCont);
	}
	myTrie->firstCont = NULL;
	if (myTrie->range) { delete [] myTrie->range; }
        myTrie->range = NULL;
}

template <class T, class S>
void destroyContainer(container<T, S>* myContainer)
{
	cout << "Destroying container" << endl;
	if (myContainer->myTrie) { destroyTrie(myContainer->myTrie); }
	myContainer->myTrie = NULL;
	if (myContainer->myTerm) { destroyTerm(myContainer->myTerm); }
	
	term<T, S> *curTerm, *tempTerm;
	curTerm = myContainer->firstTerm;
	for (int i = 0; i < myContainer->itemCount; i++)
	{
		tempTerm = curTerm;
		curTerm = curTerm->nextTerm;
		destroyTerm(tempTerm);
	}
	myContainer->firstTerm = NULL;
	myContainer->itemCount = 0;
}

template <class T, class S>
void destroyTerm(term<T, S>* myTerm)
{
	cout << "Destroying term" << endl;
	if (myTerm->coef) { delete myTerm->coef; }
	if (myTerm->exponents) { delete [] myTerm->exponents; }
	myTerm->nextTerm = NULL;
}

//creates necessary containers, if any, for insertion of myValue into myTrie
template <class T, class S>
void ensureContainers(burstTrie<T, S>* myTrie, S myValue)
{
    //cout << "Ensuring " << myValue << " is within [" << myTrie->range[0] << ", " << myTrie->range[1] << "]" << endl;
	container<T, S> *curContainer, *lastContainer;

	if (myValue < myTrie->range[0])
	{
		cout << "Setting new firstContainer" << endl;
		lastContainer = myTrie->firstCont;
		myTrie->firstCont = curContainer = createContainer<T, S>();
		
		for (S i = myValue + 1; i < myTrie->range[0]; i++)
		{
			cout << "Creating new container " << i << endl;
			curContainer->nextCont = createContainer<T, S>();
			curContainer = curContainer->nextCont;
		}
		curContainer->nextCont = lastContainer;
		myTrie->range[0] = myValue;
	}
	else if (myValue > myTrie->range[1])
	{
		cout << "Updating sortMax" << endl;
		curContainer = myTrie->firstCont;
		for (S i = myTrie->range[0]; i < myTrie->range[1]; i++)
		{ curContainer = curContainer->nextCont; }
		for (S i = myTrie->range[1] + 1; i <= myValue; i++)
		{
			cout << "Creating new container " << i << endl;
			curContainer->nextCont = createContainer<T, S>();
			curContainer = curContainer->nextCont;
		}
		myTrie->range[1] = myValue;
	}
}

//inserts a copy of myTerm into a myTrie
template <class T, class S>
void trieInsert(burstTrie<T, S>* myTrie, term<T, S>* myTerm)
{   
	cout << "Adding element to trie!" << endl;
        cout << *myTerm->coef;
        if (myTerm->expCount > 0)
        {
            cout << ", [";
            for (int i = 0; i < myTerm->expCount; i++)
            { cout << myTerm->exponents[i] << ", "; }
            cout << "]" << endl;
            ensureContainers(myTrie, myTerm->exponents[0]);
        }
        

	container<T, S>* myContainer = myTrie->firstCont;
        int j = 0;
        for (S i = myTrie->range[0]; i < myTerm->exponents[0]; i++)
	{
	    if (myContainer->nextCont == NULL)
	    { cout << "Warning: next element is null, current container " << i << endl; }
	    myContainer = myContainer->nextCont; j++;
	}
        cout << "Inserting into container " << j << endl;
        term<T, S>* newTerm;
        newTerm = stripTerm<T, S>(myTerm);
        
        cout << "Stripped term: ";
        cout << *newTerm->coef;
        if (newTerm->expCount > 0)
        {
            cout << ", [";
            for (int i = 0; i < newTerm->expCount; i++)
            { cout << newTerm->exponents[i] << ", "; }
            cout << "]" << endl;
        }
        
        if (newTerm->expCount == 0)
        {
            if(myContainer->myTerm)
	    {
		*newTerm->coef += *myContainer->myTerm->coef;
	    }
	    myContainer->myTerm = newTerm;
        }
	else if (myContainer->itemCount >= 0 && myContainer->itemCount + 1 < BURST_MAX) //regular container insert
	{
            containerInsert(myContainer, newTerm);
	}
	else //inserting into a trie instead, might have to create it from the container
	{
	    if ( myContainer->itemCount > 0 ) { containerToTrie(myContainer); }
	    container<T, S>* curContainer = myContainer->myTrie->firstCont;
	    ensureContainers(myContainer->myTrie, newTerm->exponents[0]);
	    for (S i = myContainer->myTrie->range[0]; i < newTerm->exponents[0]; i++)
	    {
		if (curContainer->nextCont == NULL)
		{ cout << "ERROR: Next element is null, trying to reach container " << i << endl; return; }
		curContainer = curContainer->nextCont;
            }	
	    trieInsert<T, S>(myContainer->myTrie, newTerm);
            destroyTerm(newTerm);
	}
	
}

//converts a container to a trie
template <class T, class S>
void containerToTrie(container<T, S>* myContainer)
{
	cout << "Converting container (" << myContainer->itemCount << " items) to trie" << endl;
	myContainer->myTrie = (burstTrie<T, S>*) malloc (sizeof(burstTrie<T, S>));
	burstTrie<T, S>* curTrie = myContainer->myTrie; //will be sorted on the first exponent of the container's terms
	
	//get min, max
	term<T, S> *nextTerm, *curTerm = myContainer->firstTerm;
        curTrie->range = new S[2];
	curTrie->range[0] = curTerm->exponents[0];
	cout << "Setting sortMin to " << curTrie->range[0];
	for (int i = 1; i < myContainer->itemCount; i++)
	{
		if (curTerm->nextTerm == NULL)
		{ cout << "Next element is null, trying to reach last term" << endl; cin; return; }
		curTerm = curTerm->nextTerm;
	}
	curTrie->range[1] = curTerm->exponents[0]; //the last element should have the highest first exponent
	cout << ", sortMax to " << curTrie->range[1] << endl;
	cout << "There are " << curTrie->range[1] - curTrie->range[0] + 1 << " containers in the new trie" << endl;
	
	//make containers
	container<T, S>* curContainer;
	curContainer = curTrie->firstCont = createContainer<T, S>();
	cout << "Created first container" << endl;
	for (S i = curTrie->range[0]; i < curTrie->range[1]; i++)
	{
		cout << "Creating new container " << i << endl;
		curContainer->nextCont = createContainer<T, S>();
		curContainer = curContainer->nextCont;
	}
	curContainer->nextCont = NULL;
	
	cout << "Containers created, inserting terms" << endl;
	//insert terms into containers
	curTerm = myContainer->firstTerm;

	for (int i = 0; i < myContainer->itemCount; i++)
	{
		cout << "Parsing term " << i << endl;
		curContainer = curTrie->firstCont;
		for (S i = curTrie->range[0]; i < curTerm->exponents[0]; i++)
		{
			if (curContainer->nextCont == NULL)
			{
				cout << "ERROR: Next element is null, trying to reach container " << i << endl;
				return;
			}
			curContainer = curContainer->nextCont;
		}		
                trieInsert(curTrie, curTerm);
		nextTerm = curTerm->nextTerm;
		destroyTerm(curTerm);
		curTerm = nextTerm;
	}
	
	myContainer->firstTerm = NULL; //using myTrie from now on
	myContainer->itemCount = -1; //not relevant for tries
}

//inserts a term into container
template <class T, class S>
void containerInsert(container<T, S>* myContainer, term<T, S>* myTerm)
{
	cout << "Inserting element into container w/ " << myContainer->itemCount << " items" << endl;
	if (myContainer->itemCount == 0)
	{
                cout << "Setting first term" << endl;
		myTerm->nextTerm = NULL;
		myContainer->firstTerm = myTerm;
		myContainer->itemCount++;
	}
	else if (myContainer->itemCount == -1) //inserting into trie
	{
		container<T, S>* curContainer = myContainer->myTrie->firstCont;
		for (S i = myContainer->myTrie->range[0]; i < myTerm->exponents[0]; i++)
		{
			if (curContainer->nextCont == NULL)
			{
				cout << "ERROR: Next element is null, trying to reach container " << i << endl;
				return;
			}
			curContainer = curContainer->nextCont;
		}
                trieInsert(myContainer->myTrie, myTerm);
		destroyTerm(myTerm);
	}
	else
	{
		int sortV;
		sortV = compareTerm(myTerm, myContainer->firstTerm);
		//find where to insert, sort on remaining variables
                if (sortV < 0)
                {
                    cout << "new firstTerm" << endl;
                    myTerm->nextTerm = myContainer->firstTerm;
                    myContainer->firstTerm = myTerm;
		    myContainer->itemCount++;
                }
                else if (sortV == 0)
		{
		    cout << "modifying firstTerm";
		    *myContainer->firstTerm->coef += (*myTerm->coef);
		    cout << " to " << *myContainer->firstTerm->coef << endl;
		}
		else
		{
                    term<T, S> *lastTerm, *curTerm;
                    lastTerm = myContainer->firstTerm; curTerm = myContainer->firstTerm->nextTerm;
                    int i = 0;
                    //while (curTerm && compareTerm(curTerm, myTerm) < 0) { lastTerm = curTerm; curTerm = curTerm->nextTerm; i++; }
                    //if (compareTerm(lastTerm, myTerm) == 0) { cout << "Equals found." << endl; }
                    //else { cout << "Inserted at position " << i << endl; }
		    while (curTerm)
		    {
			sortV = compareTerm(curTerm, myTerm);
			cout << "Compared, result is " << sortV << endl;
			if (sortV > 0) { cout << "Breaking out" << endl; break; }
			if (sortV == 0)
			{
				cout << "Combining matching terms: ";
				*curTerm->coef += *myTerm->coef;
				destroyTerm(myTerm);
				cout << *curTerm->coef << endl;
				break;
			}
			lastTerm = curTerm;
			curTerm = curTerm->nextTerm; i++;
		    }
		    if (sortV != 0)
		    {
			lastTerm->nextTerm = myTerm;
			myTerm->nextTerm = curTerm;
			myContainer->itemCount++;
		    }
                }
	}
}

//returns -1 if first < second, 1 if first > second
template <class T, class S>
int compareTerm(term<T, S>* firstTerm, term<T, S>* secondTerm)
{
        if (firstTerm->degree != -1 && secondTerm->degree != -1)
        {
            if (firstTerm->degree > secondTerm->degree) { return 1; }
            else if (firstTerm->degree < secondTerm->degree) { return -1; }
        }
	for (int i = 0; i < firstTerm->expCount && i < secondTerm->expCount; i++)
	{
		if (firstTerm->exponents[i] < secondTerm->exponents[i]) { return -1; }
		if (firstTerm->exponents[i] > secondTerm->exponents[i]) { return 1; }
	}
	if (firstTerm->expCount == secondTerm->expCount && firstTerm->degree == secondTerm->degree)
	{ cout << "Note: comparing identical arrays" << endl; return 0; }
	else if (firstTerm->expCount < secondTerm->expCount) //shorter exp array is 'smaller'
	{ return -1; }
	else
	{ return 1; }
}

template <class T, class S>
class TrieIterator {
public:
  // Take term and consume it.
  virtual void consumeTerm(term<T, S>* myTerm) = 0;
  void init () {}
  void setDimension(int myVal) { dimension = myVal; if (curPrefix) { /*delete [] curPrefix;*/ } curPrefix = new S[dimension]; curDepth = 0; }
  int getDimension() { return dimension; }
  burstTrie<T, S>* myTrie;
  int curDepth; //0 is root burst trie, 1 is its child, etc.
  S* curPrefix; //to store the last (curDepth + 1) values that are implied by burst trie ordering
  int dimension;

  term<T, S>* getFullTerm(term<T, S>* thisTerm)
  {
	term<ZZ, int>* newTerm = (term<ZZ, int>*) malloc(sizeof(term<ZZ, int>));
	newTerm->coef = new ZZ((*thisTerm->coef));
	cout << "New coef is " << *newTerm->coef << endl;
	newTerm->expCount = dimension;
	cout  << "exp count is " << newTerm->expCount << endl;
	newTerm->exponents = new int[dimension];
	cout << "Depth is " << curDepth << endl;
	int i = 0;
	for (; i < curDepth; i++) //sorted values
	{
		cout << "11. Adding " << curPrefix[i] << endl;
		newTerm->exponents[i] = curPrefix[i];
	}
	for (; i < thisTerm->expCount + curDepth; i++) //stored exps
	{
		cout << "12. Adding " << thisTerm->exponents[i - curDepth] << endl;
		newTerm->exponents[i] = thisTerm->exponents[i - curDepth];
	}
	for ( ; i < dimension; i++) //trailing zeroes
	{
		cout << "13. Adding 0" << endl;
		newTerm->exponents[i] = 0;
	}
        newTerm->degree = -1;
	return newTerm;
  }

  void enumContainer(container<T, S>* myContainer)
  {
        cout << "Enumerating container at depth " << curDepth << " (" << myContainer->itemCount << " items)" << endl;
            
        term<T, S>* curTerm = myContainer->firstTerm;
            
        for (int i = 0; i < myContainer->itemCount; i++)
        {
            consumeTerm(curTerm);
            curTerm = curTerm->nextTerm;
        }
  } 
  
  void enumTrie(burstTrie<T, S>* curTrie)
  {
        cout << "Enumerating trie at depth " << curDepth << endl;
        cout << "Range is [" << curTrie->range[0] << ", " << curTrie->range[1] << "]" << endl;
        container<T, S>* curContainer = curTrie->firstCont;
              
        for (S i = curTrie->range[0]; i <= curTrie->range[1]; i++)
        {
            curPrefix[curDepth] = i;
            curDepth++; 
            if (curContainer && curContainer->myTerm) //print coefficient in trie
            {
                consumeTerm(curContainer->myTerm);
            }      
            if (curContainer && curContainer->itemCount == -1) //this is a trie, test if myTrie is NULL instead of itemCount?
            {
                enumTrie(curContainer->myTrie);
            }
            else if (curContainer && curContainer->itemCount > 0) //regular container
            {
                enumContainer(curContainer);
            }
            curDepth--;
            cout << "Enumerated container " << i << endl;
            curContainer = curContainer->nextCont;
        }
  }
};

template <class T, class S>
class TriePrinter: public TrieIterator<T, S> {
    public:
    stringstream output;
    
    TriePrinter()
    {
        stringstream output (stringstream::in | stringstream::out);
    }
    
    void consumeTerm(term<T, S>* myTerm)
    {
        if (output.str() != "") { output << ", "; }
        if (myTerm->degree > -1)
        {
            output << "[" << *myTerm->coef << ", [" << myTerm->degree;
        }
        else
        {
            output << "[" << *myTerm->coef;
        }
        output << ", [" << TrieIterator<T, S>::curPrefix[0];
	int i = 1;
	for (; i < TrieIterator<T, S>::curDepth; i++) //sorted values
	{
		output << ", " << TrieIterator<T, S>::curPrefix[i];
	}
	for (; i < myTerm->expCount + TrieIterator<T, S>::curDepth; i++) //stored exps
	{
		output << ", " << myTerm->exponents[i - TrieIterator<T, S>::curDepth];
	}
	for ( ; i < TrieIterator<T, S>::dimension; i++) //trailing zeroes
	{
		output << ", 0";
	}
        /*for (int i = 1; i < TrieIterator<T, S>::curDepth; i++)
        {
            output << ", " << TrieIterator<T, S>::curPrefix[i];
        }

	for (int i = 0; i < myTerm->expCount; i++) //prints the stored exponents
	{
		output << ", " << myTerm->exponents[i];
	}
	for (int i = TrieIterator<T, S>::curDepth + myTerm->expCount + 1; i < TrieIterator<T, S>::dimension; i++) //prints trailing zero exponents
	{
		output << ", 0";
	}*/
	output << "]]";
        if (myTerm->degree > -1) { output << "]"; }
    }
    
    string printTrie(burstTrie<T, S>* trie)
    {
        output.str("");
        enumTrie(trie);
        return output.str();
    }
};

#endif
