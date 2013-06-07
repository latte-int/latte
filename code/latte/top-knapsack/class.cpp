

#include <iostream>
#include <string>

using namespace std;

class Star
{
public:
	string s;
	Star()
	{
	
	}
	~Star()
	{
		cout << "Good buy " << s.c_str() << endl;
	}
	void setName(const char *t)
	{
		s = string(t);
	}
	
	const Star & operator=(const Star & rhs)
	{
		if (this == &rhs)
			return *this;
		s = rhs.s;
		return *this;
	}
	
	const Star  operator*(const Star & rhs) const
	{
		Star t;
		t.setName((s+rhs.s).c_str());
		return t;
	}
};


void foo(const Star & g, const Star & h)
{
	Star t;
	t = g*h;
	cout << "This is t saying my new name is " << t.s.c_str() << endl;
	
}
int main(void)
{
	Star g,h;
	g.setName("Brandon");
	h.setName("Bryan");
	foo(g, h);


	return 0;
}
