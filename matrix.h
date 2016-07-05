/*
 * matrix.h
 * Copyright (C) 2016 zobekenobe <zobekenobe@gmail.com>
 *
 * matrix.h is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * matrix.h is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef matrix_h
#define matrix_h

#include"vector.h"

template<typename T> class matrix;

template<typename T> const matrix<T> operator+(const matrix<T>&, const matrix<T>&);
template<typename T> const matrix<T> operator-(const matrix<T>&, const matrix<T>&);
template<typename T> const matrix<T> operator*(const matrix<T>&, const matrix<T>&);
template<typename T> const matrix<T> operator,(const matrix<T>&, const matrix<T>&);

template<typename T> const matrix<T> operator+ (const matrix<T>&, const T&);
template<typename T> const matrix<T> operator+ (const T&, const matrix<T>&);
template<typename T> const matrix<T> operator- (const matrix<T>&, const T&);
template<typename T> const matrix<T> operator- (const T&, const matrix<T>&);
template<typename T> const matrix<T> operator* (const matrix<T>&, const T&);
template<typename T> const matrix<T> operator* (const T&, const matrix<T>&);
template<typename T> const matrix<T> operator/ (const matrix<T>&, const T&);


template<typename T> const vector<T> operator* ( vector<T>&,  matrix<T>&);
template<typename T> const vector<T> operator* ( matrix<T>&,  vector<T>&);
template<typename T> const matrix<T> operator^ ( vector<T>&,  vector<T>&);

template<typename T> std::istream& operator>>(std::istream&, const matrix<T>&);
template<typename T> std::ostream& operator<<(std::ostream&, const matrix<T>&);

template<typename T> matrix<T> ones(int, int);
template<typename T> matrix<T> zeros(int,int);
template<typename T> matrix<T> eye(int);
template<typename T> void swap(T&, T&);

template<typename T>
class matrix
{
	friend vector<T>;

	private:
		int dim[2];
		vector<T>** elem;
	public:
		matrix();
		matrix(int);
		matrix(int , int);
		matrix(const matrix&);
		~matrix();

		// Operator Overloading
		// Binary Operators: matrix - matrix operations

		friend const matrix<T> operator+ <T>(const matrix<T>&, const matrix<T>&);
		friend const matrix<T> operator- <T>(const matrix<T>&, const matrix<T>&);
		friend const matrix<T> operator* <T>(const matrix<T>&, const matrix<T>&);
		friend const matrix<T> operator, <T>(const matrix<T>&, const matrix<T>&);
		
		// Binary Operations: matrix - constant operations
		friend const matrix<T> operator+ <T>(const matrix<T>&, const T&);
		friend const matrix<T> operator+ <T>(const T&, const matrix<T>&);
		friend const matrix<T> operator- <T>(const matrix<T>&, const T&);
		friend const matrix<T> operator- <T>(const T&, const matrix<T>&);
		friend const matrix<T> operator* <T>(const matrix<T>&, const T&);
		friend const matrix<T> operator* <T>(const T&, const matrix<T>&);
		friend const matrix<T> operator/ <T>(const matrix<T>&, const T&);

		friend std::ostream& operator<< <T>(std::ostream&,const matrix<T>&);
		friend std::istream& operator>> <T>(std::istream&,const matrix<T>&);
		// Unary Operator
		const matrix<T>& operator =(const matrix<T>&);
		const matrix<T>& operator+=(const matrix<T>&);
		const matrix<T>& operator-=(const matrix<T>&);
		const matrix<T>& operator*=(const matrix<T>&);
		
		const matrix<T>& operator+=(const auto&);
		const matrix<T>& operator-=(const auto&);
		const matrix<T>& operator*=(const auto&);
		const matrix<T>& operator/=(const auto&);	
			
		const matrix<T>& operator++();
		const matrix<T>& operator++(int);
		
		const vector<T>  operator[](int);
		
		int* __();	
			
		T& operator()(int, int);	
		inline  int _r(){return dim[0];}
		inline  int _c(){return dim[1];}
		vector<T>& row(int i){return elem[i][0];}
		
		// Matrix Generator 
		friend matrix<T> ones <T>(int,int);
		friend matrix<T> zeros <T>(int, int);
		friend matrix<T> eye <T>(int);
		friend void swap <T>(T&,T&);

		// post-operables
		vector<T> _slicerow(int);
		vector<T> _slicecol(int);
		matrix<T> _t(); 
		
		friend const vector<T> operator* <T>( matrix<T>&, vector<T>&); 	// v = m.u
 		friend const vector<T> operator* <T>( vector<T>&, matrix<T>&);	// v = u.m	
 		friend const matrix<T> operator^ <T>( vector<T>&, vector<T>&);	// m = u.v
 		
 		
		//template<typename U>  friend T norm<T>(const vector<T>&);		// norm
		//template<typename U>	friend T norm<T>(const matrix<T>&);
		//template<typename U>	friend T norminf<T>(const vector<T>&);		// infinite norm
		//template<typename U>	friend T norminf<T>(const matrix<T>&);
};

#endif 
