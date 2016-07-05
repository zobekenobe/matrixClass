/*
 * SLE.cpp
 * Copyright (C) 2016 zobekenobe <zobekenobe@gmail.com>
 *
 * SLE.cpp is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * SLE.cpp is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include<cmath>

template<typename T>
void swapp(T& a, T& b)
{
	T temp = a;a = b;b = temp;
}

template<typename T>
matrix<T> eyes(int n)
{
	matrix<T> e(n);
	for(size_t i=0;i<n;i++)
		e(i,i) = 1.0;
	return e;
}

template<typename T>
void gauss0(matrix<T>& g,vector<T>& v)
{ 
	int n = g._r();
	for(int i=0; i<n-1;i++)
	{
		for(int j=i+1;j<n;j++)
		{
			T constant  = g(j,i)/g(i,i);
			for(int k=0;k<n;k++)
			{
				g(j,k)-= g(i,k)*constant;
			}   v[j]  -= v(i)*constant;
		}
	}
}

template<typename T>
matrix<T> pivot(matrix<T> m)
{
	//finds the maximum of a given column and pivots 
	//to have a the diagonals with the maximum valued 
	//element. returns a vector of arranged indices
	int n = m._r();
	matrix<T> p(std::move(eyes<T>(n)));
	for(int col=0;col<n;col++)
	{
		int max = col;
		for(int i=col;i<n;i++)
		{
			if(fabs(m(max,col)<m(i,col)))
				max = i;
		}
		swapp(m.row(max),m.row(col));
		swapp(p.row(max),p.row(col));
	}
	return p;
}

template<typename T>
vector<T> bsub(matrix<T>& m, vector<T>& v)
{
	int n = v._l();
	vector<T> res(n);
	res[n-1] = v(n-1)/m(n-1,n-1);
	for(int i=n-2;i>=0;i--)
	{
		T sum {};
		for(int j=n-1;j>i;j--)
		{
			sum += m(i,j)*res(j);	
		} res[i] = (v(i)-sum)/m(i,i);
	} 
	return res;
}

template<typename T>
vector<T> gaussElim0(matrix<T>& m, vector<T>& v)
{
	gauss0(m,v);
	return bsub(m,v);
}

template<typename T>
vector<T> gaussElimP(matrix<T>& m, vector<T>& v)
{
	matrix<T> p(std::move(pivot(m)));
	m = p*m;
	v = p*v;
	gauss0(m,v);
	return bsub(m,v);
}

template<typename T>
matrix<T> inv(matrix<T> m)
{
	int r = m._r();
	matrix<T> I(std::move(eyes<T>(r)));
	for(int i=0;i<r-1;i++)
	{
		for(int j=i+1;j<r;j++)
		{
			T constant = m(j,i)/m(i,i);
			for(int k=0;k<r;k++)
			{
				I(j,k) -= I(i,k)*constant;
				m(j,k) -= m(i,k)*constant;
			}
		}
	}
	for(int i=r-1;i>0;i--)
		for(int j=i-1;j>=0;j--)
		{
			T constant = m(j,i)/m(i,i);
			for(int k=0;k<r;k++)
				{
					I(j,k) -= I(i,k)*constant;
					m(j,k) -= m(i,k)*constant;
				}	
		}

	for(int i=0;i<r;i++)
		I.row(i) /= m(i,i);

	return I;
}

template<typename T>
void gaussJor(matrix<T>& m, vector<T>& v)
{
	int r = m._r();
	for(int i=0;i<r-1;i++)
	{
		for(int j=i+1;j<r;j++)
		{
			T constant = m(j,i)/m(i,i);
			for(int k=0;k<r;k++)
				m(j,k)-= m(i,k)*constant;
			v[j] -= v[i]*constant;
		}
	}

	for(int i=r-1;i>0;i--)
	{
		for(int j=i-1;j>=0;j--)
		{
			T constant = m(j,i)/m(i,i);
			for(int k=0;k<r;k++)
				m(j,k) -= m(i,k)*constant;
			v[j] -= v[i]*constant;
		}
	}

	for(int i=0;i<r;i++)
		v[i]/=m(i,i);
}
/*	LU Decomposition (Factorization): Triangulation
*	LU Decomposition (facorization) of a non-singular (square) matrix A
* 	means expressing the matrix as the multiplication of a lower triangular
*	matrix L and an upper triangular matrix U, where a lower/upper triangular
*	matrix is a matrix having no nonzero elements above/below the diagonal. For
*	the case where row switching operation is needed like in the Gauss 
*	elimination, we include a permutation matrix P representing the necessary 
*	row switching operations to write the LU decomposition as PA = LU */

template<typename T>
matrix<T> LU_crout(matrix<T>& m)
{
	matrix<T> p(std::move(pivot(m)));
	//m = p*m;
	//v = p*v;
	matrix<T> lu(m._r());
	int n = m._r();
	for(int i=0;i<n;i++)
	{
		lu(i,0) = m(i,0);
		for(int j=1;j<n;j++)
		{
			if(j<=i)
			{	
				T sum {};
				for(int k=0;k<j;k++)
				{
					sum += lu(i,k)*lu(k,j);
				} lu(i,j) = (m(i,j)-sum);
			}
			else
			{	
				T sum {};
				for(int k=0;k<i;k++)
				{
					sum += lu(i,k)*lu(k,j);
				} lu(i,j) = ((m(i,j)-sum)/lu(i,i));
			}
		}
	}	
	return lu;
}

template<typename T>
vector<T> fsub(matrix<T>& m, vector<T>& v)
{
	vector<T> ans(v._l());
	for(int i=0;i<m._r();i++)
	{
		T sum {};
		for(int j=0;j<i;j++)
		{
			sum+= m(i,j)*ans[j];
		}ans[i] = (v[i]-sum)/m(i,i); 
	}
	return ans;
}

template<typename T>
vector<T> crout(matrix<T> m, vector<T>& v)
{
	/*
	*	The Crout method simplifies the whole factorization by setting
	* 	the values of the diagonal elements of matrix U equal to 1, 
	*	which has the advantage over the Gauss Elimination method where
	* 	the number of unknowns in L and U is n^2-n, which is lower
	*/
	vector<T> ans(v._l());
	matrix<T> lu(m._r());
	lu = LU_crout(m);
	ans = fsub(lu,v);
	for(int i=0;i<lu._r();i++)
		lu(i,i) = 1.0;
	ans = bsub(lu,ans);
	return ans;
}

template<typename T>
vector<T> doolittle(matrix<T> m, vector<T> v)
{ 
	/*
	*	The Doolittle method is the other side of the Crout, where the 
	*	values of the disgonal elements of the lower triangular matrix
	*	is set to 1. 
	*/
	vector<T> ans(v._l());
	matrix<T> lu(m._r());
	lu = LU_crout(m); lu = lu._t();
	for(int i=0;i<lu._r();i++)
		lu(i,i) = 1.0;
	ans = bsub(lu,ans);
	ans = fsub(lu,v);

	return ans;
}

template<typename T>
matrix<T> LU_cholesky(matrix<T> m)
{
	matrix<T> ch(m._r());
	for(int i=0;i<m._r();i++)
	{
		for(int j=0;j<=i;j++)
		{
			T sum {};
			if(i==j)
			{
				for(int k=0;k<j;k++)
					sum+=ch(j,k)*ch(j,k);
				ch(j,j) = sqrt(m(j,j)-sum);
			}
			else
			{
				for(int k=0;k<j;k++)
					sum+= ch(i,k)*ch(j,k);
				ch(i,j) = (m(i,j)-sum)/ch(j,j);
				//std::cout<<m(j,j)<<" ";
				//std::cout<<ch(i,j)<<std::endl;
			}
		}
	}
	return ch;
}


template<typename T>
vector<T> cholesky(matrix<T> m, vector<T> v)
{
	vector<T> ans(v._l());
	matrix<T> lu(m._r());
	lu = LU_cholesky(m);
	ans = fsub(lu,v);
	lu = lu._t();
	ans = bsub(lu,ans);
	return ans;
}

template<typename T>
matrix<T> LU_tridag(matrix<T>& m)
{
	matrix<T> lu(m._r());
	int n = m._r()-1;
	lu(0,0) = m(0,0);lu(0,1) = m(0,1)/lu(0,0);
	for(int i=1;i<(m._r()-1);i++)
	{
		lu(i,i-1) = m(i,i-1);
		lu(i,i) = m(i,i) - (lu(i,i-1)*lu(i-1,i));
		lu(i,i+1) = m(i,i+1)/lu(i,i);
	}
	lu(n,n-1) = m(n,n-1);
	lu(n,n) = m(n,n)-(lu(n,n-1)*lu(n-1,n));

	return lu;
}

template<typename T>
vector<T> tridag(matrix<T> m, vector<T> v)
{
	matrix<T> lu(m._r());
	lu = LU_tridag(m);
	vector<T> ans(v._l());
	ans = fsub(lu,v);
	for(int i=0;i<lu._r();i++)
		lu(i,i) = 1.0;	
	ans = bsub(lu,ans);
	return ans;
}

template<typename T>
matrix<T> band(matrix<T>& m, int m1, int m2)
{
	matrix<T> b(m._r(),m1+m2+1);
	for(int i=0;i<m._r();i++)
	{
		b(i,m1) = m(i,i);
		for(int j=0;j<m1;j++)
		{
			b(i,j) = ((i+j)<m1 )?0:m(i,i-(m1-j));
		}
		for(int j=1;j<=m2;j++)
		{
			b(i,m1+j) = ((i+j)>=m._r())?0:m(i,i+j);
		}
	}
	return b;
}