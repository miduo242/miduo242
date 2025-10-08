/*
 *   This file is part of the OpenPhase (R) software library.
 *
 *   Copyright (c) 2009-2022 Ruhr-Universitaet Bochum,
 *                 Universitaetsstrasse 150, D-44801 Bochum, Germany
 *             AND 2018-2022 OpenPhase Solutions GmbH,
 *                 Universitaetsstrasse 136, D-44799 Bochum, Germany.
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *   File created :   2011
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev;
 *                         Philipp Engels; Raphael Schiedung
 *
 */


#ifndef MATRIX_H
#define MATRIX_H

#include <cassert>

#include "Base/Includes.h"

namespace openphase
{

template <class T>
class OP_EXPORTS Matrix /// Matrix template class. Can handle any type of values
{
 public:
    Matrix()
    {
        Array = nullptr;
        Size_X = 0;
        Size_Y = 0;
    }
    Matrix(const Matrix<T>& rhs)
    {
        Size_X = rhs.size_X();
        Size_Y = rhs.size_Y();

        Array = new T[Size_X*Size_Y] ();

        for (size_t i = 0; i < Size_X; i++)
        for (size_t j = 0; j < Size_Y; j++)
        {
            Array[Size_Y*i+j] = rhs(i,j);
        }
    };

    Matrix(const size_t nx, const size_t ny)
    {
        Allocate(nx, ny);
    }
    T& operator()(const size_t x, const size_t y) const
    {
        assert(x < Size_X && "Access beyond storage range");
        assert(y < Size_Y && "Access beyond storage range");

        return Array[Size_Y*x + y];
    }
    void Allocate(const size_t nx, const size_t ny)
    {
        Size_X = nx;
        Size_Y = ny;
        Array = new T[nx*ny] ();
    }
    void Reallocate(const size_t nx, const size_t ny)
    {
        if (Array != nullptr) delete[] Array;
        Size_X = nx;
        Size_Y = ny;
        Array = new T[nx*ny] ();
    }
    void set(const size_t x, const size_t y, const T value)
    {
        assert(x < Size_X && "Access beyond storage range");
        assert(y < Size_Y && "Access beyond storage range");

        Array[Size_Y*x + y] = value;
    }
    void set_to_value(const T value)
    {
        for (size_t x = 0; x < Size_X; x++)
        for (size_t y = 0; y < Size_Y; y++)
        {
            Array[Size_Y*x + y] = value;
        }
    }
    void set_to_zero()
    {
        for (size_t x = 0; x < Size_X; x++)
        for (size_t y = 0; y < Size_Y; y++)
        {
            Array[Size_Y*x + y] = T();
        }
    }
    void add(const size_t x, const size_t y, const T value)
    {
        assert(x < Size_X && "Access beyond storage range");
        assert(y < Size_Y && "Access beyond storage range");

        Array[Size_Y*x + y] += value;
    }
    T& get(const size_t x, const size_t y) const
    {
        assert(x < Size_X && "Access beyond storage range");
        assert(y < Size_Y && "Access beyond storage range");

        return Array[Size_Y*x + y];
    }
    size_t size_X() const
    {
        return Size_X;
    }
    size_t size_Y() const
    {
        return Size_Y;
    }
    T get_min() const
    {
        T min = std::numeric_limits<T>::max();
        for (size_t x = 0; x < Size_X; x++)
        for (size_t y = 0; y < Size_Y; y++)
        {
            if (Array[Size_Y*x + y] < min) min = Array[Size_Y*x + y];
        }
        return min;
    }
    T get_max() const
    {
        T max = std::numeric_limits<T>::min();
        for (size_t x = 0; x < Size_X; x++)
        for (size_t y = 0; y < Size_Y; y++)
        {
            if (Array[Size_Y*x + y] > max) max = Array[Size_Y*x + y];
        }
        return max;
    }
    bool IsNotAllocated() const
    {
        return (Array == nullptr);
    }
    bool IsAllocated() const
    {
        return !(Array == nullptr);
    }
    std::string print(void) const
    {
        std::stringstream out;
        for(size_t i = 0; i < Size_X; i++)
        for(size_t j = 0; j < Size_Y; j++)
        {
            out << "||" << std::setprecision(6) << std::right
                        << std::setw(8) << Array[Size_Y*i+j];
            if (j == Size_Y-1)
            {
                out << "||\n";
            }
            else
            {
                out << " ";
            }
        }
        return out.str();
    };

    Matrix<T>& invert(void)
    {
        Matrix<T> tmp = this->inverted();
        for(size_t i = 0; i < Size_X; i++)
        for(size_t j = 0; j < Size_Y; j++)
        Array[Size_Y*i+j] = tmp.get(i,j);
        return *this;
    };

    bool is_invertable(void)
    {
        if (Size_X != Size_Y)
        {
            //Matrix is not square
            return false;
        }
        std::vector<std::vector<T> > Out;
        std::vector<std::vector<T> > Uni;
        std::vector<T> temp;
        temp.resize(Size_X, 0);
        for(size_t i = 0; i < Size_X; i++)
        {
            Out.push_back(temp);
            Uni.push_back(temp);
        }

        std::vector<size_t> indxc(Size_X,0);
        std::vector<size_t> indxr(Size_X,0);
        std::vector<size_t> ipiv(Size_X,0);
        size_t icol = 0;
        size_t irow = 0;
        T pivinv;
        T dum;
        for(size_t i = 0; i < Size_X; i++)
        {
            for(size_t j = 0; j < Size_Y; j++)
            {
                if (i == j)
                {
                    Uni[i][j] = 1.0;
                }
                else
                {
                    Uni[i][j] = 0.0;
                }
                Out[i][j] = Array[Size_Y*i+j];
            }
        }

        for(size_t i = 0; i < Size_X; i++)
        {
            T big = 0.0;
            for(size_t j = 0; j < Size_X; j++)
            if(ipiv[j] != 1)
            for(size_t k = 0; k < Size_Y; k++)
            {
                if(ipiv[k] == 0)
                {
                    if(fabs(Out[j][k]) >= big)
                    {
                        big = fabs(Out[j][k]);
                        irow = j;
                        icol = k;
                    };
                }
                else if (ipiv[k] > 1)
                {
                    return false;
                }
            };
            ++(ipiv[icol]);
            if(irow != icol)
            {
                for (size_t l = 0; l < Size_X; l++)
                {
                    T temp = Out[irow][l];
                    Out[irow][l] = Out[icol][l];
                    Out[icol][l] = temp;
                };
                for (size_t l = 0; l < Size_X; l++)
                {
                    T temp = Uni[irow][l];
                    Uni[irow][l] = Uni[icol][l];
                    Uni[icol][l] = temp;
                };
            };
            indxr[i] = irow;
            indxc[i] = icol;

            if (fabs(Out[icol][icol]) <= 0.0/*DBL_EPSILON*/)
            {
                return false;
            }
            pivinv = 1.0/Out[icol][icol];
            Out[icol][icol] = 1.0;
            for(size_t l = 0; l < Size_X; l++) Out[icol][l] *= pivinv;
            for(size_t l = 0; l < Size_X; l++) Uni[icol][l] *= pivinv;
            for(size_t ll = 0; ll < Size_X; ll++)
            if(ll != icol)
            {
                dum = Out[ll][icol];
                Out[ll][icol] = 0.0;
                for(size_t l = 0; l < Size_X; l++) Out[ll][l] -= Out[icol][l]*dum;
                for(size_t l = 0; l < Size_X; l++) Uni[ll][l] -= Uni[icol][l]*dum;
            }
        }
        return true;
    };

    Matrix<T> inverted(void)
    {
        if (Size_X != Size_Y)
        {
            std::cerr << "Matrix: Can Not Compute Inverse Matrix.\n"
                      << this->print() << "Matrix is not quadratic!!!\n";
            std::terminate();
        }
        std::vector<std::vector<T> > Out;
        std::vector<std::vector<T> > Uni;
        std::vector<T> temp;
        temp.resize(Size_X, 0);
        for(size_t i = 0; i < Size_X; i++)
        {
            Out.push_back(temp);
            Uni.push_back(temp);
        }

        std::vector<size_t> indxc(Size_X,0);
        std::vector<size_t> indxr(Size_X,0);
        std::vector<size_t> ipiv(Size_X,0);
        size_t icol = 0;
        size_t irow = 0;
        T pivinv;
        T dum;
        for(size_t i = 0; i < Size_X; i++)
        {
            for(size_t j = 0; j < Size_Y; j++)
            {
                if (i == j)
                {
                    Uni[i][j] = 1.0;
                }
                else
                {
                    Uni[i][j] = 0.0;
                }
                Out[i][j] = Array[Size_Y*i+j];
            }
        }

        for(size_t i = 0; i < Size_X; i++)
        {
            T big = 0.0;
            for(size_t j = 0; j < Size_X; j++)
            if(ipiv[j] != 1)
            for(size_t k = 0; k < Size_Y; k++)
            {
                if(ipiv[k] == 0)
                {
                    if(fabs(Out[j][k]) >= big)
                    {
                        big = fabs(Out[j][k]);
                        irow = j;
                        icol = k;
                    };
                }
                else if (ipiv[k] > 1)
                {
                    std::cerr << "dMatrix6x6: Can Not Compute Inverse Matrix.\n"
                              << this->print() << "Matrix is Singular 1!!!\n";
                    std::terminate();
                }
            };
            ++(ipiv[icol]);
            if(irow != icol)
            {
                for (size_t l = 0; l < Size_X; l++)
                {
                    T temp = Out[irow][l];
                    Out[irow][l] = Out[icol][l];
                    Out[icol][l] = temp;
                };
                for (size_t l = 0; l < Size_X; l++)
                {
                    T temp = Uni[irow][l];
                    Uni[irow][l] = Uni[icol][l];
                    Uni[icol][l] = temp;
                };
            };
            indxr[i] = irow;
            indxc[i] = icol;

            if (fabs(Out[icol][icol]) <= DBL_EPSILON)
            {
                std::cerr << "Matrix: Can Not Compute Inverse Matrix.\n"
                          << this->print() << "Matrix is Singular 2!!!\n";
                std::terminate();
            }
            pivinv = 1.0/Out[icol][icol];
            Out[icol][icol] = 1.0;
            for(size_t l = 0; l < Size_X; l++) Out[icol][l] *= pivinv;
            for(size_t l = 0; l < Size_X; l++) Uni[icol][l] *= pivinv;
            for(size_t ll = 0; ll < Size_X; ll++)
            if(ll != icol)
            {
                dum = Out[ll][icol];
                Out[ll][icol] = 0.0;
                for(size_t l = 0; l < Size_X; l++) Out[ll][l] -= Out[icol][l]*dum;
                for(size_t l = 0; l < Size_X; l++) Uni[ll][l] -= Uni[icol][l]*dum;
            }
        }
        for(int l = Size_X-1; l >= 0; l--)
        {
            if(indxr[l] != indxc[l])
            for(size_t k = 0; k < Size_X; k++)
            {
                T temp = Out[k][indxr[l]];
                Out[k][indxr[l]] = Out[k][indxc[l]];
                Out[k][indxc[l]] = temp;
            };
        }

        Matrix<T> temp2;
        temp2.Allocate(Size_X, Size_Y);
        for(size_t i = 0; i < Size_X; i++)
        for(size_t j = 0; j < Size_Y; j++)
        {
             temp2.set(i,j,Out[i][j]);
        }
        return temp2;
    };

    /*Matrix<T>& transpose(void)
    {
        Matrix<T> tmp = transposed();

        Reallocate(tmp.size_X(),tmp.size_Y());
        for(size_t i = 0; i < Size_X; i++)
        for(size_t j = 0; j < Size_Y; j++)
        Array[Size_Y*i+j] = tmp.get(i,j);
        return *this;
    };*/

    Matrix<T> transposed(void)
    {
        Matrix<T> temp;
        temp.Allocate(Size_Y,Size_X);

        for(size_t i = 0; i < temp.size_X(); i++)
        for(size_t j = 0; j < temp.size_Y(); j++)
        {
            temp.set(j,i,Array[temp.size_Y()*i+j]);
        }

        return temp;
    };

    Matrix<T> operator*(const Matrix<T>& rhs) const
    {
        assert(Size_X == rhs.size_X() && "Matrices do not have the same dimension!");
        assert(Size_Y == rhs.size_Y() && "Matrices do not have the same dimension!");

        size_t size_X = rhs.size_Y();
        size_t size_Y = Size_X;

        std::vector<T> tmp(size_Y, 0.0);
        std::vector<std::vector<T> > Out(size_X,tmp);

        for (size_t i = 0; i < size_X; i++)
        for (size_t j = 0; j < size_Y; j++)
        {
            double value = 0.0;
            for (size_t r = 0; r < size_Y; r++)
            {
                value += get(i,r)*rhs.get(r,j);
            }
            Out[i][j] = value;
        }

        Matrix<T> tmp2(size_X,size_Y);

        for (size_t i = 0; i < size_X; i++)
        for (size_t j = 0; j < size_Y; j++)
        {
            tmp2.Array[size_Y*i+j] = Out[i][j];
        }

        return tmp2;
    };

    Matrix<T> operator+(const Matrix<T>& rhs)
    {
        assert(Size_X == rhs.size_X() && "Matrices do not have the same dimension!");
        assert(Size_Y == rhs.size_Y() && "Matrices do not have the same dimension!");

        Matrix<T> Out = *this;

        for (size_t i = 0; i < Size_X; i++)
        for (size_t j = 0; j < Size_Y; j++)
        {
            Out(i,j) += rhs(i,j);
        }
        return Out;
    };

    Matrix<T>& operator+=(const Matrix<T>& rhs)
    {
        assert(Size_X == rhs.size_X() && "Matrices do not have the same dimension!");
        assert(Size_Y == rhs.size_Y() && "Matrices do not have the same dimension!");

        for (size_t i = 0; i < Size_X; i++)
        for (size_t j = 0; j < Size_Y; j++)
        {
            Array[Size_Y*i+j] += rhs(i,j);
        }
        return *this;
    };

    Matrix<T>& operator=(const Matrix<T>& rhs)
    {
        assert(Size_X == rhs.size_X() && "Matrices do not have the same dimension!");
        assert(Size_Y == rhs.size_Y() && "Matrices do not have the same dimension!");

        for (size_t i = 0; i < Size_X; i++)
        for (size_t j = 0; j < Size_Y; j++)
        {
            Array[Size_Y*i+j] = rhs(i,j);
        }
        return *this;
    };

    ~Matrix()
    {
        delete[] Array;
    }
    T* data()
    {
        return Array;
    }

    static Matrix<T> max(const Matrix<T>& M1, const Matrix<T>& M2)
    {
        assert(M1.size_X() == M2.size_X() && "Matrices do not have the same dimension!");
        assert(M1.size_Y() == M2.size_Y() && "Matrices do not have the same dimension!");

        Matrix<T> Out(M1.size_X(), M1.size_Y());
        for(size_t n = 0; n < M1.size_X(); n++)
        for(size_t m = 0; m < M1.size_Y(); m++)
        {
            Out(n,m) = std::max(M1(n,m), M2(n,m));
        }
        return Out;
    }

    static Matrix<T> min(const Matrix<T>& M1, const  Matrix<T>& M2)
    {
        assert(M1.size_X() == M2.size_X() && "Matrices do not have the same dimension!");
        assert(M1.size_Y() == M2.size_Y() && "Matrices do not have the same dimension!");

        Matrix<T> Out(M1.size_X(), M1.size_Y());
        for(size_t n = 0; n < M1.size_X(); n++)
        for(size_t m = 0; m < M1.size_Y(); m++)
        {
            Out(n,m) = std::min(M1(n,m), M2(n,m));
        }
        return Out;
    }

 protected:
 private:
    T*     Array;
    size_t Size_X;
    size_t Size_Y;
};

}// namespace openphase
#endif
