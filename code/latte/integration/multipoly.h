/*This code is to store multivariate polynomials by total degree. The structure head contains a linked list that 
stores all monomials with the same degree. The structure head itself uses a linked list. Node stores a monomial with integer
coefficient and at most 50 degrees. nil and nill are used as dummy variables to make our code simpler.*/


#include <stdlib.h>
#include <iostream>
using namespace std;

short compare(const short x[50],const short y[50], short top)
//compare two numerical strings of the same degree by lexicographic order
{              
        for (short i=0;i<top;i++)
        {
            if (x[i]>y[i]) {return 1;};
            if (x[i]<y[i]) {return -1;};
        };
        return 0;
};       
class multipoly
{
      public:
             short num_variable;
             struct node{        //a node represents a monomial
                    short expo[50];     
                    long coef;
                    node* next;
                        };
             struct head{                    //head is a collection of monomials with the same total degree
                    short degree;
                    head* link;
                    node* nill;              //nill is a dummy variable
                    };
             head* nil;
             
             
      multipoly(string inputpoly)
      {               //constructor , which takes in a string of prescribed format
             nil=new head;
             nil->link=nil;
             nil->nill=new node;
             (nil->nill)->next=nil->nill;
             short index=1;
             short l=inputpoly.length();
             num_variable=-1;
             while (index<l-1) 
             {
                  node* x=new node;
                  int loc_comma=inputpoly.find(",",index);
                  string temp=inputpoly.substr(index+1,loc_comma-1-index);
                  x->coef=atol(temp.c_str());
                  short deg=0;
                  int loc_bra=inputpoly.find("]",loc_comma);
                  temp=inputpoly.substr(loc_comma+2,loc_bra-2-loc_comma);
                  if (temp.length()==0) {num_variable=0;}
                  else 
                  {
                   int t=0;
                   int tt=temp.find(",",0);
                   int counter=0;
                   string subtemp;
                   while (tt!=string::npos) 
                   {
                        subtemp=temp.substr(t,tt-t);
                        counter++;
                        x->expo[counter-1]=atoi(subtemp.c_str());
                        deg+=x->expo[counter-1];
                        t=tt+1;tt=temp.find(",",t);
                   }
                   subtemp=temp.substr(t,temp.length()-t);
                   counter++;
                   x->expo[counter-1]=atoi(subtemp.c_str());
                   deg+=x->expo[counter-1];
                   if (num_variable==-1) num_variable=counter;
                  }
                  index=loc_bra+3;
                  head *temp1=nil;                                                        //start insertion  
                  head *temp2=nil->link;                                                  
                  while ((temp2!=nil)&&(deg>temp2->degree))
                  {
                        temp1=temp2;
                        temp2=temp2->link;
                  }
                  if ((temp2==nil)||(deg<temp2->degree)) {                                //create a new collection
                                  head *h=new head; 
                                  temp1->link=h;
                                  h->link=temp2;
                                  h->degree=deg;
                                  h->nill=new node;
                                  (h->nill)->next=x;
                                  x->next=h->nill;
                                  }
                  else {                                                                   //insert a monomial into temp2
                         node *p=temp2->nill;
                         node *q=p->next;
                         short t=compare(x->expo,q->expo,num_variable);
                         while ((q!=temp2->nill)&&(t==1))
                                        {
                                            p=q;
                                            q=q->next;
                                            if (q!=temp2->nill) {t=compare(x->expo,q->expo,num_variable);};
                                        }
                         if ((q==temp2->nill)||(t<0))                                       //insert a new monomial between p,q
                         {
                             p->next=x;
                             x->next=q;
                         }
                         else 
                         {
                            q->coef+=x->coef;
                            if (q->coef==0) {
                                            p->next=q->next;free(q);
                                            if ((temp2->nill)->next==temp2->nill) {temp1->link=temp2->link;free(temp2);};
                                            }
                         };
                       };
             }
      }
      ~multipoly(){}                                        
      void print()                              //printing polynomial by convention
             {
                  head *tt=nil->link;
                  if (tt==nil) {cout<<"0"<<endl; return;}; 
                  while (tt!=nil) 
                  { 
                       bool init=1;
                       cout<<"degree"<<tt->degree<<":";
                       node *ttt=(tt->nill)->next;
                       while (ttt!=tt->nill)
                       {
                             if (!init) 
                             {
                                if (ttt->coef>0) cout<<"+";
                                if (ttt->coef==-1) {cout<<"-";} else 
                                if (ttt->coef!=1) {cout<<ttt->coef;};
                             }
                             else 
                             {
                                  if (ttt->coef==-1) {cout<<"-";} else
                                  if (ttt->coef!=1) {cout<<ttt->coef;};
                                  init=0;
                             };
                             for (int j=0;j<num_variable;j++) {
                             if (ttt->expo[j]!=0) {cout<<"x"<<j+1; if (ttt->expo[j]>1) cout<<"^"<<ttt->expo[j];};
                                                };    
                             if ((tt->degree==0)&&((ttt->coef==1)||(ttt->coef==-1))) {cout<<1;};  
                             ttt=ttt->next;
                       }  
                       tt=tt->link;cout<<endl;
                  }
             }              
};
