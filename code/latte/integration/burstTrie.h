#ifndef BURSTTRIE_H
#define BURSTTRIE_H

#define BURST_MAX 10
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <stdio.h>
#include <sstream>
#include <assert.h>

#define BT_DEBUG 1

NTL_CLIENT

struct trieElem
{
    bool isTrie;
    void* myVal;
    
    trieElem* next;
};

template <class T, class S>
class BurstTerm
{
public:
    BurstTerm(int myLength)
    {
        length = myLength;
	exps = NULL;
    }
    
    BurstTerm(const T& newCoef, S* newExps, int start, int myLength, int myDegree)
    {
	//cout << "Creating term: ";
        degree = myDegree;
        length = myLength - start;
        exps = new S[length];
        for (int i = start; i < myLength; i++)
	{ exps[i - start] = newExps[i]; }
        coef = newCoef;
        next = NULL;
    }
    
    ~BurstTerm()
    {
	//cout << "Destroying term" << endl;
        if (length > 0 && exps)
        { delete [] exps; }
    }
    
    bool lessThan(BurstTerm<T, S>* other, bool &equal)
    {
        equal = false;
        for (int i = 0; i < length && i < other->length; i++)
        {
	    //cout << "Comparing " << exps[i] << " v. " << other->exps[i] << endl;
            if (exps[i] < other->exps[i]) { return true; }
            if (exps[i] > other->exps[i]) { return false; }
        }
        if (length < other->length) { return true; }
        if (length > other->length) { return false; }
        equal = true;
        return false;
    }
    
    BurstTerm<T, S>* next;
    
    T coef;
    S* exps;
    int length;
    int degree;
};

template <class T, class S> class BurstTrie;

template <class T, class S>
class BurstCont
{
public:
    BurstCont()
    {
	//cout << "New container" << endl;
        firstTerm = NULL;
        termCount = 0;
    }
    
    ~BurstCont()
    {
	//cout << "Destroying container (" << termCount << " terms) ..." << endl;
        BurstTerm<T, S> *temp, *old;
        temp = firstTerm;
        for (int i = 0; i < termCount; i++)
        {
            old = temp->next;
            delete temp;
            temp = old;
        }
    }
    
    void insertTerm(const T& newCoef, S* newExps, int start, int myLength, int myDegree)
    {
	//cout << "Inserting term into container" << endl;
        if (firstTerm == NULL)
        {
            firstTerm = new BurstTerm<T, S>(newCoef, newExps, start, myLength, myDegree);
            termCount++;
            return;
        }
        
	bool equal;
        BurstTerm<T, S>* newTerm = new BurstTerm<T, S>(newCoef, newExps, start, myLength, myDegree);
        if (newTerm->lessThan(firstTerm, equal))
        {
            newTerm->next = firstTerm;
            firstTerm = newTerm;
            termCount++;
            return;
        }
	if (equal)
        {
            firstTerm->coef += newTerm->coef;
            delete newTerm;
            return;
        }
        
        BurstTerm<T, S> *curTerm = firstTerm;
        BurstTerm<T, S> *oldTerm;
        
        while(curTerm && curTerm->lessThan(newTerm, equal))
        {
            oldTerm = curTerm;
            curTerm = curTerm->next;
        }
	if (equal)
        {
            curTerm->coef += newTerm->coef;
            delete newTerm;
            return;
        }
        
        if (curTerm == NULL) 
        {
            oldTerm->next = newTerm;
        }
        else //oldTerm < newTerm < curTerm
        {
            oldTerm->next = newTerm;
            newTerm->next = curTerm;
        }
        termCount++;        
    }
    
    BurstTrie<T, S>* burst()
    {
        BurstTrie<T, S>* myTrie = new BurstTrie<T, S>();
        BurstTerm<T, S>* curTerm = firstTerm;
        BurstTerm<T, S>* oldTerm;
        for (int i = 0; i < termCount; i++)
        {
            myTrie->insertTerm(curTerm->coef, curTerm->exps, 0, curTerm->length, curTerm->degree);
            oldTerm = curTerm->next;
            curTerm = oldTerm;
        }
        return myTrie;
    }
    
    BurstTerm<T, S>* getTerm(int index)
    {
        assert(index < termCount);
        BurstTerm<T, S>* myTerm = firstTerm;
        for (int i = 0; i < index; i++)
        {  myTerm = myTerm->next; }
        return myTerm;
    }
    
    int termCount;
private:
    BurstTerm<T, S>* firstTerm;
};

template <class T, class S>
class BurstTrie
{
public:
    BurstTrie()
    {
        curIndex = -1;
        range = NULL;
        firstElem = NULL;
    }
    
    ~BurstTrie()
    {
	//cout << "Destroying trie" << endl;
        if (range)
        {
            trieElem *temp = firstElem;
            trieElem *old;
	     while (temp != NULL)
            {
		//cout << "Destroying trie element.." << endl;
                //destroy element container or trie
                if (temp->isTrie)
                {
                    delete (BurstTrie<T, S>*)temp->myVal;
                }
                else
                {
                    delete (BurstCont<T, S>*)temp->myVal;
                }
                old = temp->next;
                //destroy trieElem
                free(temp); 
                temp = old;
            }
            delete [] range;
        }
    }
    
    void insertTerm(const T& newCoef, S* newExps, int start, int myLength, int myDegree)
    {
        assert(myLength > 0);
	/*cout << "Inserting term into trie: " << newCoef;
	for (int i = start; i < myLength; i++)
	{
	    cout << ", " << newExps[i];
	}
	cout << endl;*/
        if (range == NULL)
        {
            range = new S[2];
            range[0] = range[1] = newExps[0];
            firstElem = (trieElem*)malloc(sizeof(trieElem));
            firstElem->next = NULL;
            firstElem->myVal = new BurstCont<T, S>();
            firstElem->isTrie = false;
        }
        else
        {
            checkRange(newExps[start]);
        }
        
        trieElem *curElem = firstElem;
        for (S i = range[0]; i < newExps[start]; i++)
        { curElem = curElem->next; }
        
        if (curElem->isTrie)
        {
            ((BurstTrie<T, S>*)curElem->myVal)->insertTerm(newCoef, newExps, start + 1, myLength, myDegree);
        }
        else
        {
	    BurstCont<T, S>* temp = (BurstCont<T, S>*)curElem->myVal;
	    //cout << "Trie element is a container (" << temp->termCount << " elements)" << endl;
            if (temp->termCount == BURST_MAX && myLength > 1)
            {
		//cout << "Bursting container..." << endl;
                BurstTrie<T, S>* newTrie = temp->burst();
		//cout << "Burst trie created, deleting container" << endl;
                delete temp;
    /*BurstTerm<ZZ, int>* abcd = new BurstTerm<ZZ, int>(myLength - start - 1);
    newTrie->begin();
    int i = 0;
    while (newTrie->nextTerm(abcd))
    {
        cout << "Term " << i++ << ": " << abcd->coef;
        for (int j = 0; j < abcd->length; j++)
	{
		cout << ", " << abcd->exps[j];
	}
	cout << endl;
    }
    delete abcd;*/
                curElem->isTrie = true;
                curElem->myVal = newTrie;
                newTrie->insertTerm(newCoef, newExps, start + 1, myLength, myDegree);
            }
            else
            {
                temp->insertTerm(newCoef, newExps, start + 1,  myLength, myDegree);
            }
        }
    }
    
    void checkRange(const S& myVal)
    {
        if (myVal < range[0]) //new minimum
        {
	    trieElem *temp = (trieElem*)malloc(sizeof(trieElem)); //new first element for myVal
	    trieElem *old = temp;
            temp->next = NULL;
            temp->myVal = new BurstCont<T, S>();
	    temp->isTrie = false;
            for (S i = myVal + 1; i < range[0]; i++)
            {
                //create new element
                temp->next = (trieElem*)malloc(sizeof(trieElem));
                //advance to it
                temp = temp->next;
                //set pointer
                temp->next = NULL;
                //allocate container
                temp->myVal = new BurstCont<T, S>();
                temp->isTrie = false;
            }
            temp->next = firstElem;
            //set new first element
            firstElem = old;
            range[0] = myVal;
        }
        else if (myVal > range[1]) //new maximum
        {
            trieElem *temp = firstElem;
            for (S i = range[0]; i < range[1]; i++)
            { temp = temp->next; }
            for (S i = range[1]; i < myVal; i++)
            {
                temp->next = (trieElem*)malloc(sizeof(trieElem));
                temp = temp->next;
                temp->next = NULL;
                temp->myVal = new BurstCont<T, S>();
                temp->isTrie = false;
            }
            range[1] = myVal;
        }
    }
    
    void begin()
    {
        curIndex = -1;
	trieElem *temp = firstElem;
        while (temp != NULL)
        {
            if (temp->isTrie)
            {
                ((BurstTrie<T, S>*)temp->myVal)->begin();
            }
            temp = temp->next;
        }
    }
    
    bool nextTerm(BurstTerm<T, S>* newTerm)
    {
        if (curIndex == -1)
        { curIndex = range[0]; termIndex = 0; }
        
        if (curIndex > range[1])
        {
            return false;
        }
        
        S* myExps = new S[newTerm->length];
        BurstTerm<T, S>* storedTerm = getTerm(myExps, 0); 
        if (storedTerm == NULL)
	{
	    //cout << "Advancing container (max is " << range[1] << ")" << endl;
	    delete[] myExps;
	    termIndex = 0;
	    do
	    {
		curIndex++;
		//cout << "Checking container " << curIndex << endl;
	    }
	    while (curIndex <= range[1] && !nextTerm(newTerm));
	    if (curIndex > range[1])
	    { return false; }
	    else
	    { return true; }
	}
        newTerm->coef = storedTerm->coef;
	newTerm->degree = storedTerm->degree;      
        
        for (int i = 0; i < storedTerm->length; i++)
        {
	    myExps[newTerm->length - storedTerm->length + i] = storedTerm->exps[i];
        }
	if (newTerm->exps) { delete [] newTerm->exps; }
	newTerm->exps = myExps;
        newTerm->degree = storedTerm->degree;
        return true;
    }
    
    BurstTerm<T, S>* getTerm(S* myExps, int myDepth)
    {
        trieElem *temp = firstElem;
	BurstTerm<T, S>* storedTerm;
        for (S i = range[0]; i < curIndex; i++)
        { temp = temp->next; }
        
	//cout << "Checking trie element" << endl;
        if (temp->isTrie)
        {
            //update prefix
	    //cout << "Updating prefix to " << curIndex << endl;
	    myExps[myDepth] = curIndex;
	    //cout << "Fetching stored term at index " << ((BurstTrie<T, S>*)temp->myVal)->curIndex << endl;
	    //check iterators
	    if (((BurstTrie<T, S>*)temp->myVal)->curIndex == -1)
	    {
		((BurstTrie<T, S>*)temp->myVal)->curIndex = ((BurstTrie<T, S>*)temp->myVal)->range[0];
		//cout << "Resetting child trie to " << ((BurstTrie<T, S>*)temp->myVal)->curIndex << endl;
		((BurstTrie<T, S>*)temp->myVal)->termIndex = 0;
	    }
	    //get stored term
            storedTerm = ((BurstTrie<T, S>*)temp->myVal)->getTerm(myExps, myDepth + 1);
	    if (storedTerm == NULL)
	    {
		curIndex++; termIndex = 0;
		if (curIndex <= range[1])
		{ return getTerm(myExps, myDepth); }
		else
		{ return NULL; }
	    }
            return storedTerm;
        }
        else
        {
	    //cout << "Found container at trie element " << curIndex << endl;
            myExps[myDepth] = curIndex;
            if (termIndex < ((BurstCont<T, S>*)temp->myVal)->termCount)
            {
		//cout << "Fetching stored term " << termIndex << endl;
		storedTerm = ((BurstCont<T, S>*)temp->myVal)->getTerm(termIndex);
		termIndex++; return storedTerm;
	    }
            else
            {
		curIndex++; termIndex = 0;
		if (curIndex <= range[1])
		{ return getTerm(myExps, myDepth); }
		else
		{ return NULL; }
	    }
        }
    }

private:
    S* range; //S can be a class or a primitve
    trieElem *firstElem; //first element in the trie
    
    //iteration markers
    int termIndex;
    S curIndex;
};

#endif
