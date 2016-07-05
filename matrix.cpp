/*
 * matrix.cpp
 * Copyright (C) 2016 zobekenobe <zobekenobe@gmail.com>
 *
 * matrix.cpp is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * matrix.cpp is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

template<typename T>
matrix<T>::matrix(){}

template<typename T>
matrix<T>::~matrix()
{
	for(int i=0;i<dim[0];i++)
		delete elem[i];
	delete [] elem;
}


template<typename T>
matrix<T>::matrix(int a):dim{a,a},elem(new vector<T>*[a])
{
	for(size_t i=0;i<a;i++)
		elem[i] = new vector<T>(a);
}

template<typename T>
matrix<T>::matrix(int a, int b):dim{a,b},elem(new vector<T>*[a])
{
	for(size_t i=0;i<a;i++)
		elem[i] = new vector<T>(b);	
}

template<typename T>
matrix<T>::matrix(const matrix<T>& m):elem(new vector<T>*[m.dim[0]])
{
	memcpy(dim,m.dim,sizeof dim);
	for(size_t i=0;i<m.dim[0];i++)
	{	elem[i] = new vector<T>(m.dim[1]);
 		elem[i][0] = m.elem[i][0];
	}
}

template<typename T>
T& matrix<T>::operator()(int a, int b)
{
	return elem[a][0][b];			//return (*(*(elem+a)+0))[b];   Alternative way of writing it
}

// Binary Operators: matrix - matrix operations

template<typename T>
const matrix<T> operator+(const matrix<T>& m1, const matrix<T>& m2)
{
	matrix<T> m(m1.dim[0],m1.dim[1]);
	for(size_t i=0;i<m.dim[0];i++)
		m.elem[i][0] = m1.elem[i][0]+m2.elem[i][0];
	return m;
}

template<typename T>
const matrix<T> operator-(const matrix<T>& m1, const matrix<T>& m2)
{
	matrix<T> m(m1.dim[0],m1.dim[1]);
	for(size_t i=0;i<m.dim[0];i++)
		m.elem[i][0] = m1.elem[i][0]-m2.elem[i][0];
	return m;
}


template<typename T>
const matrix<T> operator*(const matrix<T>& m1, const matrix<T>& m2)
{
	matrix<T> m(m1.dim[0],m1.dim[1]);
	for(int i=0;i<m.dim[0];i++){
		for(int j=0;j<m.dim[1];j++){
			T sum {};
			for(int k=0;k<m.dim[0];k++){
				sum += m1.elem[i][0][k]*m2.elem[k][0][j];}
			m.elem[i][0][j] = sum;}}
	return m;
}

template<typename T>
const matrix<T> operator,(const matrix<T>& m1, const matrix<T>& m2)
{
	matrix<T> m(m1.dim[0],m1.dim[1]);
	for(int i=0;i<m.dim[0];i++)
		 m.elem[i][0] = (m1.elem[i][0],m2.elem[i][0]);
	return m;
}


// Binary Operators: Matrix - constant operations
template<typename T>
const matrix<T> operator+(const matrix<T>& m1, const T& k )
{
	matrix<T> m(m1.dim[0],m1.dim[1]);
	for(int i=0;i<m1.dim[0];i++)
		m.elem[i][0] = m1.elem[i][0] + k;
	return m;
} 

template<typename T>
const matrix<T> operator+(const T& k, const matrix<T>& m1)
{
	return (m1+k);
}

template<typename T>
const matrix<T> operator-(const matrix<T>& m1, const T& k)
{
	matrix<T> m(m1.dim[0],m1.dim[1]);	
	for(int i=0;i<m1.dim[0];i++)
		m.elem[i][0] = m1.elem[i][0] - k;
	return m;
}

template<typename T>
const matrix<T> operator-(const T& k, const matrix<T>& m1)
{
	matrix<T> m(m1.dim[0],m1.dim[1]);
	for(int i=0;i<m1.dim[0];i++)
		m.elem[i][0] = k-m1.elem[i][0];
	return m;
}

template<typename T>
const matrix<T> operator*(const matrix<T>& m1, const T& k)
{
	matrix<T> m(m1.dim[0],m1.dim[1]);
	for(int i=0;i<m1.dim[0];i++)
		m.elem[i][0] = m1.elem[i][0] * k;
	return m;
}

template<typename T>
const matrix<T> operator*(const T& k, const matrix<T>& m1)
{
	return (m1*k);
}

template<typename T>
const matrix<T> operator/(const matrix<T>& m1, const T& k )
{
	matrix<T> m(m1.dim[0],m1.dim[1]);
	for(int i=0;i<m1.dim[0];i++)
		m.elem[i][0] = m1.elem[i][0] / k;
	return m;
}
// Unary Operators:

template<typename T>
const matrix<T>& matrix<T>::operator=(const matrix<T>& m)
{
	memcpy(dim,m.dim,sizeof dim);
	for(int i=0;i<dim[0];i++)
		elem[i][0] = m.elem[i][0];
	return *this;
}

template<typename T>
const matrix<T>& matrix<T>::operator+=(const matrix<T>& m1)
{
	for(int i=0;i<dim[0];i++)
		elem[i][0] += m1.elem[i][0];
	return *this;
}

template<typename T>
const matrix<T>& matrix<T>::operator-=(const matrix<T>& m1)
{
	for(int i=0;i<dim[0];i++)
		elem[i][0] -= m1.elem[i][0];
	return *this;
}

template<typename T>
const matrix<T>& matrix<T>::operator*=(const matrix<T>& m1)
{
	for(int i=0;i<dim[0];i++)
		elem[i][0] *= m1.elem[i][0];
	return *this;
}

template<typename T>
const matrix<T>& matrix<T>::operator+=(const auto& k)
{
	for(int i=0;i<dim[0];i++)
		elem[i][0] += k;
	return *this;
}

template<typename T>
const matrix<T>& matrix<T>::operator-=(const auto& k)
{
	for(int i=0;i<dim[0];i++)
		elem[i][0] -= k;
	return *this;
}

template<typename T>
const matrix<T>& matrix<T>::operator*=(const auto& k)
{
	for(int i=0;i<dim[0];i++)
		elem[i][0] *= k;
	return *this;
}

template<typename T>
const matrix<T>& matrix<T>::operator/=(const auto& k)
{
	for(int i=0;i<dim[0];i++)
		elem[i][0] /= k;
	return *this;
}

template<typename T>
const matrix<T>& matrix<T>::operator++()
{
	for(int i=0;i<dim[0];i++)
		++elem[i][0];
	return *this;
}

template<typename T>
const matrix<T>& matrix<T>::operator++(int)
{
	for(int i=0;i<dim[0];i++)
		elem[i][0]++;
	return *this;
}

template<typename T>
const vector<T> matrix<T>::operator[](int i)
{
	return elem[i][0];
}

// Matrix Vector operations
template<typename T>
const vector<T> operator*( matrix<T>& m, vector<T>& v)
{
	vector<T> u(v._l());
	for(int i=0;i<v._l();i++)
	{
		u[i] = m.elem[i][0]*v;
	}
	return u;
}

template<typename T>
const vector<T> operator*( vector<T>& v, matrix<T>& m)
{
	m = m._t(); 
	return (m*v);
}

template<typename T>
const matrix<T> operator^(vector<T>& v, vector<T>& u)
{
	matrix<T> m(v._l(),u._l());
	for(size_t i=0;i<m.dim[0];i++)
		m.elem[i][0] = v(i)*u;
	return m;
}


// matrix generators
template<typename T>
matrix<T> ones(int a, int b)
{
	return (matrix<T>(a,b)+1.0);
}

template<typename T>
matrix<T> eye(int a)
{
	matrix<T> m(a);
	for(size_t i=0;i<m.dim[0];i++)
		m(i,i) = 1.0;
}

template<typename T>
void swap(T& a, T& b)
{
	T temp = a;
	a = b;
	b = temp;	
}

template<typename T>
matrix<T> matrix<T>::_t()
{
	matrix<T> m(dim[1],dim[0]);
	for(size_t i=0;i<dim[1];i++)
		for(size_t j=0;j<dim[0];j++)
			m(i,j) = this->operator()(j,i);
	return m;
}


// read and write functions
template<typename T>
std::istream& operator>>(std::istream& in, const matrix<T>& m)
{
	for(size_t i=0;i<m.dim[0];i++)
		in>>m.elem[i][0];   
	return in;
}


template<typename T>
std::ostream& operator<<(std::ostream& out,const matrix<T>& m)
{ 
	for(size_t i=0;i<m.dim[0];i++)
		out<<m.elem[i][0]<<std::endl;   //out<<*(m.elem[i]);		Alternative way of writing it
	return out;
}




