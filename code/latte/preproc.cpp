#include <iostream.h>
#include <stdio.h>
#include <fstream>
#include <NTL/ZZ.h>
#include <NTL/config.h>
#include <string.h>
#include <gmp.h>
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>

#define MAX_DIM 1000
#define HugeInt ZZ
#define LEQ '~'
#define GEQ '`'

void reset_coefficients(HugeInt *c, int max_subscript);  
void preprocess(ifstream &in, int &max_subscript, int &num_lines);

void syntax_error(int current_line);

void write_first_line(ofstream &out, int max_subscript, int num_lines);
void write_relation(ofstream &out, char relation, HugeInt *c, int max_subscript);
void write_last_line(ofstream &out, int *equation_lines, int equ_count);

int get_coeff(ifstream &in, HugeInt &c);
int get_subscript(ifstream &in, int &s);
int get_symbol(ifstream &in, char &s);

int is_relation(char s);
int char_to_digit(char d);

void skip_white(ifstream &in);


int main(int argc, char *argv[])
{
	//char *buffer = new char[10000];
	HugeInt *coefficients = new HugeInt[MAX_DIM + 1];
	//int *signs = new int[MAX_DIM + 1];
	int equ_count = 0, current_line = 0;
	int *equation_lines = new int[10000];

	int subscript, max_subscript, num_lines;
	HugeInt coeff;
	char relation, symbol;

	if(argc < 3)
	{
		cout << "\nERROR:  Too few arguments!  Proper usage is  " << argv[0] << " sourcefile targetfile\n";
		return 1;
	}

	ifstream in_pre(argv[1]);
	ofstream out(argv[2]);

	if(in_pre.fail() || out.fail())
	{
		cout << "\nFile Error!\n";
		return 2;
	}

	reset_coefficients(coefficients, MAX_DIM);

	preprocess(in_pre, max_subscript, num_lines);

	if(max_subscript == 0 || num_lines == 0)
	{
		syntax_error(current_line);
		cout << "0";
	}

	write_first_line(out, max_subscript, num_lines);

	if(max_subscript > MAX_DIM)
	{
		cout << "\nDimension too large\n";
		return 0;
	}

	in_pre.close();
	ifstream in(argv[1]);

	while(!in.eof())
	{	
		current_line++;
		relation = '0';
		symbol = '+';

		while(symbol != ';')
		{
			if(!get_coeff(in, coeff))
			{
				syntax_error(current_line);
				cout << "1";
			}

			if(!get_subscript(in, subscript))
			{
				syntax_error(current_line);
				cout << "2";
			}
					
			if(relation == '0' && symbol == '+')
				coefficients[subscript] -= coeff;
			else if( relation != '0' && symbol == '-')
				coefficients[subscript] -= coeff;
			else
				coefficients[subscript] += coeff;

			if(!get_symbol(in, symbol))
			{
				syntax_error(current_line);
				cout << "3";
			}

			if(is_relation(symbol))
			{
				if(relation != '0')
				{
					syntax_error(current_line);
					cout << "4";
				}

				else
				{
					relation = symbol;
					symbol = '+';
					if(relation == '=')
						equation_lines[equ_count++] = current_line;
				}
			}
				
		}

		if(relation == '0')
		{
			syntax_error(current_line);
			cout << "5";
		}
			
		write_relation(out, relation, coefficients, max_subscript);
		reset_coefficients(coefficients, max_subscript); 			

	}	

	write_last_line(out, equation_lines, equ_count);

	return 0;
}

void skip_white(ifstream &in)
{
	char read;

	while(1)
	{
		if(in.eof())
			break;
		read = in.peek();
		if(read != ' ' && read != '\n')
			if(read != '\t' && read != '\r')
				break;

		in.get(read);

	}

	return;
}

void reset_coefficients(HugeInt *c, int max_subscript)
{
	for(int i = 0; i <= max_subscript; i++)
		c[i] = 0;

	return;

}

/* void preprocess(ifstream &in, int &max_subscript, int &num_lines)
{
	char read;
	max_subscript = 0;
	num_lines = 0;
	int d, subscript;

	while(!in.eof())
	{
		in.get(read);

		if(read == ';')
			num_lines++; */
void preprocess(ifstream &in, int &max_subscript, int &num_lines)
{

  char read = 'x';
  max_subscript = 0;
  num_lines = 0;
  int d, subscript;

  int flag = 0;

  while(!in.eof())
    {

      in.get(read);

      if((read == ';') && (flag == 0))
	{
	  num_lines++;
	  flag = 1;
	}

      else if(read != ';')
	flag = 0;

      if(read == '[') //everything starting from this line is exactly the same
	//		if(read == '[')
		{
			subscript = 0;

			while(1)
			{
				in.get(read);
				d = char_to_digit(read);
				if(d == -1)
					break;
				subscript *= 10;
				subscript += d;
				if(in.eof())
					break;
			}

			if(subscript > max_subscript)
				max_subscript = subscript;
		}
		
	}		
				 	
	return;
}

void syntax_error(int current_line)
{
	cout << "\nSyntax error before semi-colon #" << current_line << "\n";
	return;
	exit(1);
} 

void write_first_line(ofstream &out, int max_subscript, int num_lines)
{
	out << num_lines << " " << max_subscript + 1;
	return;

}

void write_relation(ofstream &out, char relation, HugeInt *c, int max_subscript)
{
	HugeInt temp;
	
	if(relation == '>' || relation == GEQ)
		for(int i = 0; i <= max_subscript; i++)
			c[i] *= -1;

	if(relation == '>' || relation == '<')
		c[0] -= 1;

	out << "\n";
	temp = c[0];
	out << temp;

	for(int i = 1; i <= max_subscript; i++)
	{
		out << " ";
		temp = c[i];
	      	out << temp;
	}
		
	return;
}

void write_last_line(ofstream &out, int *equation_lines, int equ_count)
{
	if(equ_count == 0)
		return;

	out << "\nlinearity " << equ_count;

	for(int i = 0; i < equ_count; i++)
		out << " " << equation_lines[i];

	return;
}

int get_coeff(ifstream &in, HugeInt &c)
{
	if(in.eof())
		return 0;

	char read;
	int d, good = 0, sign = 1;

	c = 0;

	skip_white(in);

	if(in.eof())
		return 0;

	read = in.peek();

	if(read == '-')
	{
		sign = -1;
		in.get(read);
	}

	while(1)
	{
		skip_white(in);

		if(in.eof())
			return 0;

		read = in.peek();

		d = char_to_digit(read);
		
		if(d == -1)
		{
			if(!good)
			{
				if(read != 'x' && read != 'X')
					return 0;
				c = sign;
			}
			else
				c *= sign;

			return 1;
		}
		else
			good = 1;

		c *= 10;
		c += d;

		in.get(read);
	}

	c *= sign;
	return good;
	
}

int get_subscript(ifstream &in, int &s)
{
	
	if(in.eof())
		return 0;
	skip_white(in);
	if(in.eof())
		return 0;

	
	char read;
	int d, good = 0;

	s = 0;



	read = in.peek();
		
	if(read != 'x' && read != 'X')
		return 1;

	in.get(read);
		


	read = in.peek();
		
	if(read != '[')
		return 0;

	in.get(read);
		
	
	while(1)
	{
		skip_white(in);

		if(in.eof())
			return 0;

		read = in.peek();
		
		d = char_to_digit(read);
		
		if(d == -1)
		{
			if(read == ']')
			{
				in.get(read);
				skip_white(in);
				return good;
			}
			else
				return 0;
		}
		else
			good = 1;

		s *= 10;
		s += d;
	 
		in.get(read);
	}

	return good;

}

int get_symbol(ifstream &in, char &s)
{
	if(in.eof())
		return 0;

	char read, read2;

	skip_white(in);
	in.get(read);
	skip_white(in);
	
	if(read == '>' || read == '<')
	{
		read2 = in.peek();

		if(read2 == '=')
		{
			if(read == '<')
				s = LEQ;
			else
				s = GEQ;
			in.get(read2);
		}
		
		else
			s = read;
	}
	else
		s = read;

	if(s != '+' & s != '-')
	if(s != '>' & s != '<')
	if(s != GEQ & s != LEQ)
	if(s != '=' & s != ';')
		return 0;

	return 1;
}

int is_relation(char s)
{
	if(s == '<' || s == '>')
		return 1;

	if(s == LEQ || s == GEQ)
		return 1;

	if(s == '=')
		return 1;

	return 0;

}



int char_to_digit(char d)
{

	if(d == '0')
		return 0;
	else if(d == '1')
		return 1;
	else if(d == '2')
		return 2;
	else if(d == '3')
		return 3;
	else if(d == '4')
		return 4;
	else if(d == '5')
		return 5;
	else if(d == '6')
		return 6;
	else if(d == '7')
		return 7;
	else if(d == '8')
		return 8;
	else if(d == '9')
		return 9;
	else
		return -1;

}
