template <class T, class S>
class BTrieIterator : public PolyIterator<T, S>
{
public:
	BTrieIterator()
	{

	}
	
	void setTrie(BurstTrie<T, S>* trie, int myDim)
	{
		assert (myDim > 0);
		myTrie = trie;
		dimension = myDim;
		triePath = new trieElem*[dimension];
		curTerm.exps = new S[dimension];
		curTerm.length = dimension;
	}
	
	void begin()
	{
		//cout << "Starting iteration" << endl;
		curDepth = -1;
		storedTerm = NULL;
	}
	
	BurstContainer<T, S>* nextContainer()
	{
		//cout << "Next container at depth " << curDepth << endl;
		trieElem* nextElem;
		if (curDepth < 0)
		{
			curDepth++;
			nextElem = triePath[0] = myTrie->firstElem;
			curTerm.exps[0] = myTrie->range[0];
			//cout << "!exp 0: " << curTerm.exps[0] << endl;
		}
		else
		{
			nextElem = triePath[curDepth]->next;
			curTerm.exps[curDepth]++;
			while (nextElem)
			{
				//cout << "exp " << curDepth << ": " << curTerm.exps[curDepth] << endl;
				if (nextElem->isTrie)
				{
					
					//cout << "trie found" << endl;
					break;
				}
				if (((BurstContainer<T, S>*)nextElem->myVal)->termCount > 0)
				{
					//cout << "container found - " << ((BurstContainer<T, S>*)nextElem->myVal)->termCount << endl;
					break;
				}
				//cout << "Skipping trie element " << curTerm.exps[curDepth] << endl;
				nextElem = nextElem->next;
				curTerm.exps[curDepth]++;
			}
			triePath[curDepth] = nextElem;
		}
		
		if (!nextElem) //end of trie, move back up
		{
			//cout << "end of trie at depth " << curDepth << endl;
			if (curDepth == 0) { return NULL; }
			curDepth--;
			return nextContainer();
		}
		
		return firstContainer(nextElem);
	}
	
	BurstContainer<T, S>* firstContainer(trieElem* myElem)
	{
		//cout << "Looking for container at depth " << curDepth << endl;
		if (myElem->isTrie)
		{
			curDepth++;
			triePath[curDepth] = ((BurstTrie<T, S>*)myElem->myVal)->firstElem;
			curTerm.exps[curDepth] = ((BurstTrie<T, S>*)myElem->myVal)->range[0];
			return firstContainer( ((BurstTrie<T, S>*)myElem->myVal)->firstElem );
		}
		else
		{
			//cout << "Container found." << endl;
			return ((BurstContainer<T, S>*)myElem->myVal);
		}
	}
	
	term<T, S>* nextTerm()
	{
		//cout << "Next term at depth " << curDepth << endl;
		if (!storedTerm) //end of container
		{
			//cout << "Advancing container" << endl;
			BurstContainer<T, S>* curContainer = nextContainer();
			if (curContainer)
			{ storedTerm = curContainer->firstTerm; }
			else
			{ return NULL; }
			
		}

		for (int i = curDepth + 1; i < dimension; i++)
		{
			curTerm.exps[i] = storedTerm->exps[i - curDepth - 1];
		}
		curTerm.coef = storedTerm->coef;
		curTerm.degree = storedTerm->degree;
		storedTerm = storedTerm->next;
		//cout << "got term w/coef " << curTerm->coef << endl;
		return &curTerm;
	}

	term<T, S>* getTerm()
	{
		return &curTerm;
	}
	
	~BTrieIterator()
	{
		delete [] triePath;
		delete [] curTerm.exps;
	}
	
private:
	BurstTrie<T, S>* myTrie; //trie to iterate over
	term<T, S> curTerm; //shared buffer to store values
	int dimension;
	
	BurstTerm<T, S>* storedTerm; //pointer to next stored term in current container
	trieElem** triePath;
	int curDepth;
};

template <class T>
class MBlockIterator : public PolyIterator<T, int>
{
public:
	MBlockIterator()
	{
		blockIndex = 0; coeffHead =  NULL; expHead =  NULL;
	}
	
	void setLists(eBlock* eHead, cBlock<T>* cHead, int myDim, int numTerms)
	{
		assert (myDim > 0);
		dimension = myDim;
		termCount = numTerms;
		
		coeffHead = cHead; expHead = eHead;
		curTerm.exps = new int[dimension];
		curTerm.length = dimension;
		curTerm.degree = -1;
	}
	
	void begin()
	{
		blockIndex = termIndex = 0;
		curCoeff = coeffHead; curExp = expHead;
	}
	
	term<T, int>* nextTerm()
	{
		if (!curCoeff || !curExp || termIndex == termCount) { return NULL; }
		
		if (blockIndex < BLOCK_SIZE)
		{
			curTerm.coef = /*to_ZZ*/(curCoeff->data[blockIndex]);
                        for (int i = 0; i < dimension; i++)
                        {
                            curTerm.exps[i] = curExp->data[i + dimension*blockIndex];
                        }
			blockIndex++;
			termIndex++;
			return &curTerm;
		}
		else
		{
			curCoeff = curCoeff->next;
			curExp = curExp->next;
			blockIndex = 0;
			return nextTerm();
		}
	}

	term<T, int>* getTerm()
	{
		return &curTerm;
	}
	
	~MBlockIterator()
	{
		delete [] curTerm.exps;
	}
	
private:
	term<T, int> curTerm; //shared buffer to store values
	int dimension;
	int termCount;
	
	eBlock* curExp;
        cBlock<T>* curCoeff;
	
	eBlock* expHead;
        cBlock<T>* coeffHead;
	int termIndex;
	int blockIndex;
};

template <class T>
class LBlockIterator : public PolyIterator<T, ZZ>
{
public:
	LBlockIterator()
	{
		blockIndex = 0; coeffHead =  NULL; linHead = NULL;
	}
        
        void setLists(lBlock* lHead, cBlock<T>* cHead, int myDim, int numTerms)
	{
		assert (myDim > 0);
		dimension = myDim;
		termCount = numTerms;
		
		coeffHead = cHead; linHead = lHead;
		curTerm.exps = new ZZ[dimension];
		curTerm.length = dimension;
		curTerm.degree = -1;
	}
	
	void begin()
	{
		blockIndex = termIndex = 0;
		curCoeff = coeffHead; curLin = linHead;
	}
	
	term<T, ZZ>* nextTerm()
	{
		if (!curCoeff || !curLin || termIndex == termCount) { return NULL; }
		
		if (blockIndex < BLOCK_SIZE)
		{
			curTerm.coef = /* to _ Z Z */ (curCoeff->data[blockIndex]);
                        vec_ZZ myCoeffs = curLin->data[blockIndex];
                        for (int i = 0; i < dimension; i++)
                        {
                            curTerm.exps[i] = myCoeffs[i];
                        }
                        curTerm.degree = curLin->degree[blockIndex];
			blockIndex++;
			termIndex++;
			return &curTerm;
		}
		else
		{
			curCoeff = curCoeff->next;
			curLin = curLin->next;
			blockIndex = 0;
			return nextTerm();
		}
	}

	term<T, ZZ>* getTerm()
	{
		return &curTerm;
	}
	
	~LBlockIterator()
	{
		delete [] curTerm.exps;
	}
	
private:
	term<T, ZZ> curTerm; //shared buffer to store values
	int dimension;
	int termCount;
	
	lBlock* curLin;
        cBlock<T>* curCoeff;
	
	lBlock* linHead;
        cBlock<T>* coeffHead;
	int termIndex;
	int blockIndex;
};
