#ifndef FILE_AUTODIFF_h
#define FILE_AUTODIFF_h

#include <iostream>
#include <math.h>
#include <assert.h>

using namespace std;

namespace ASC_ode
{
    template <typename SCAL = double>
    class AutoDiff
    {   
        size_t size;
        SCAL val;
        VectorXd dval;
        /* SCAL ddval[D?D*D:1]={0}; */
    public:

        typedef AutoDiff<SCAL> TELEM;
        typedef SCAL TSCAL;

        /// copy constructor
        AutoDiff  (const AutoDiff & ad2) throw()
        {   
            size = ad2.size;
            val = ad2.val;
            
            for (size_t i = 0; i < size; i++) {
                dval(i) = ad2.dval(i);
            }
        }

        /// initial object with constant value
        AutoDiff  (size_t asize, SCAL aval) throw()
        {
            val = aval;
            size = asize;
            dval.resize(size);
            dval.setZero();
            for (size_t i = 0; i < size; i++)
            dval(i) = 0;
            
        }

        /// init object with (val, e_diffindex)
        AutoDiff  (size_t asize, SCAL aval, size_t diffindex)  throw()
        {
            val = aval;
            size = asize;
            dval.resize(size);
            dval.setZero();
            for (size_t i = 0; i < size; i++) {
                dval(i) = 0;
            }
            dval(diffindex) = 1;
        }

        AutoDiff (SCAL aval, const VectorXd grad)
        {
            val = aval;
            dval.resize(grad.size());
            for (size_t i = 0; i <size; i++)
            {
                dval(i) = grad(i);
            }
        }

        /// assign constant value
        AutoDiff & operator= (SCAL aval) throw()
        {
            val = aval;
            for (size_t i = 0; i < size; i++)
                dval(i) = 0;
            return *this;
        }

        /// returns value
        SCAL Value() const throw() { return val; }

        size_t Size() const throw() { return size; }

        /// returns partial derivative
        SCAL DValue (size_t i) const throw() { return dval(i); }

        VectorXd & DValue () throw() { return dval;}
 
        /// access value
        SCAL & Value() throw() { return val; }

        /// accesses partial derivative 
        SCAL & DValue (size_t i) throw() { return dval(i); }
        

        /// add AutoDiff object
        AutoDiff<SCAL> & operator+= (const AutoDiff<SCAL> & y) throw()
        {   
            assert(this->size == y.size);
            val += y.val;
            for (size_t i = 0; i < y.size; i++)
                dval(i) += y.dval(i);
            /* for (size_t i = 0; i < D*D; i++)
                ddval(i) += y.ddval(i); */
            return *this;
        }

        /// subtract AutoDiff object
        AutoDiff<SCAL> & operator-= (const AutoDiff<SCAL> & y) throw()
        {   
            assert(this->size == y.size);

            val -= y.val;
            for (size_t i = 0; i < y.size; i++)
                dval(i) -= y.dval(i);
            return *this;
        }

        /// multiply with AutoDiff object
        AutoDiff< SCAL> & operator*= (const AutoDiff< SCAL> & y) throw()
        {
            assert(this->size == y.size);

            for (size_t i = 0; i < y.size; i++)
            {
                dval(i) *= y.val;
                dval(i) += val * y.dval(i);
            }
            val *= y.val;
            return *this;
        }

        /// multiply with scalar
        AutoDiff<SCAL> & operator*= (const SCAL & y) throw()
        {
            assert(this->size == y.Size());

            for (size_t i = 0; i < y.Size(); i++) {
                dval(i) *= y;
            }
            val *= y;
            return *this;
        }

        /// divide by scalar
        AutoDiff<SCAL> & operator/= (const SCAL & y) throw()
        {
            SCAL iy = 1.0 / y;
            
            for (size_t i = 0; i < this->size; i++) {
                dval(i) *= iy;
            }
            val *= iy;
            return *this;
        }
    };

    template<typename SCAL>
    inline ostream & operator<< (ostream & ost, const AutoDiff<SCAL> & x)
    {
        ost << x.Value() << ", D = ";
        for (int i = 0; i < x.Size(); i++)
            ost << x.DValue(i) << " ";
        return ost;
    }

    template< typename SCAL, typename SCAL2, typename std::enable_if<std::is_convertible<SCAL2,SCAL>::value, int>::type = 0>
    inline AutoDiff<SCAL> operator+ (SCAL2 x, const AutoDiff<SCAL> & y) throw()
    {   
        AutoDiff<SCAL> res(y.Size(), 0);
        res.Value() = x+y.Value();
        for (size_t i = 0; i < y.Size(); i++) {
            res.DValue(i) = y.DValue(i);
        }
        /* for (size_t i = 0; i < D*D; i++) {
            res.DDValue(i) = y.DDValue(i);
        } */
        return res;
    }

    template< typename SCAL>
    inline AutoDiff< SCAL> operator+ (const AutoDiff<SCAL> & x, const AutoDiff<SCAL> & y) throw()
    {
        assert(x.Size() == y.Size());
        AutoDiff<SCAL> res(y.Size(), 0);
        res.Value() = x.Value()+y.Value();
        for (size_t i = 0; i < y.Size(); i++) {
            res.DValue(i) = x.DValue(i) + y.DValue(i);
        }

        return res;
    }

    ///
    template< typename SCAL, typename SCAL2, typename std::enable_if<std::is_convertible<SCAL2,SCAL>::value, int>::type = 0>
    inline AutoDiff< SCAL> operator+ (const AutoDiff< SCAL> & y, SCAL2 x) throw()
    {
        AutoDiff< SCAL> res(y.Size(), 0);
        res.Value() = x+y.Value();
        for (size_t i = 0; i < y.Size(); i++) {
            res.DValue(i) = y.DValue(i);
        }
        
        return res;
    }

    /*
    /// minus AutoDiff
    template< typename SCAL>
    inline AutoDiff<SCAL> operator- (const AutoDiff<SCAL> & x) throw()
    {
        AutoDiff<SCAL> res(x.Size(), 0);
        res.Value() = -x.Value();
        for (int i = 0; i < D; i++)
            res.DValue(i) = -x.DValue(i);
        return res;
    }
    */

    template< typename SCAL>
    inline AutoDiff< SCAL> operator- (const AutoDiff< SCAL> & x, const AutoDiff< SCAL> & y) throw()
    {
        assert(x.Size() == y.Size());
        AutoDiff< SCAL> res(y.Size(), 0);
        res.Value() = x.Value()-y.Value();
        for (size_t i = 0; i < y.Size(); i++) {
            res.DValue(i) = x.DValue(i) - y.DValue(i);
        }
        
        return res;
    }


    ///
    template< typename SCAL, typename SCAL2, typename std::enable_if<std::is_convertible<SCAL2,SCAL>::value, int>::type = 0>
    inline AutoDiff<SCAL> operator- (const AutoDiff<SCAL> & x) throw()
    {
        AutoDiff< SCAL> res(x.Size(), 0);
        res.Value() = -x.Value();
        for (size_t i = 0; i < x.Size(); i++) {
            res.DValue(i) = -x.DValue(i);
        }
        
        return res;
    }

    ///
    template< typename SCAL, typename SCAL2, typename std::enable_if<std::is_convertible<SCAL2,SCAL>::value, int>::type = 0>
    inline AutoDiff< SCAL> operator- (const AutoDiff< SCAL> & x, SCAL2 y) throw()
    {
        AutoDiff< SCAL> res(x.Size(), 0);
        res.Value() = x.Value()-y;
        for (size_t i = 0; i < x.Size(); i++) {
            res.DValue(i) = x.DValue(i);
        }
        
        return res;
    }

    ///
    template< typename SCAL, typename SCAL2, typename std::enable_if<std::is_convertible<SCAL2,SCAL>::value, int>::type = 0>
    inline AutoDiff< SCAL> operator- (SCAL2 x, const AutoDiff< SCAL> & y) throw()
    {
        AutoDiff< SCAL> res(y.Size(), 0);
        res.Value() = x-y.Value();
        for (size_t i = 0; i < y.Size(); i++) {
            res.DValue(i) = -y.DValue(i);
        }
        
        return res;
    }


    ///
    template< typename SCAL, typename SCAL2, typename std::enable_if<std::is_convertible<SCAL2,SCAL>::value, int>::type = 0>
    inline AutoDiff< SCAL> operator* (const SCAL2 x, const AutoDiff< SCAL> & y) throw()
    {
        AutoDiff< SCAL> res(y.Size(), 0);
        res.Value() = x*y.Value();
        for (size_t i = 0; i < y.Size(); i++) {
            res.DValue(i) = x*y.DValue(i);
        }
        
        return res;
    }

    ///
    template< typename SCAL, typename SCAL2, typename std::enable_if<std::is_convertible<SCAL2,SCAL>::value, int>::type = 0>
    inline AutoDiff< SCAL> operator* (const AutoDiff< SCAL> & y, SCAL2 x) throw()
    {
        AutoDiff< SCAL> res(y.Size(), 0);
        res.Value() = x*y.Value();
        for (size_t i = 0; i < y.Size(); i++) {
            res.DValue(i) = x*y.DValue(i);
        }
        
        return res;
    }

    ///
    template< typename SCAL>
    inline AutoDiff< SCAL> operator* (const AutoDiff< SCAL> & x, const AutoDiff< SCAL> & y) throw()
    {
        assert(x.Size() == y.Size());
        AutoDiff< SCAL> res(y.Size(), 0);
        SCAL hx = x.Value();
        SCAL hy = y.Value();

        res.Value() = hx*hy;
        for (size_t i = 0; i < x.Size(); i++) {
            res.DValue(i) = hx*y.DValue(i) + hy*x.DValue(i);
        }

        return res;
    }



    template< typename SCAL>
    inline AutoDiff< SCAL> Inv (const AutoDiff< SCAL> & x)
    {
        AutoDiff< SCAL> res(x.Size(), 1.0 / x.Value());
        for (size_t i = 0; i < x.Size(); i++) {
            res.DValue(i) = -x.DValue(i) / (x.Value() * x.Value());
        }

        return res;
    }


    template< typename SCAL>
    inline AutoDiff< SCAL> operator/ (const AutoDiff< SCAL> & x, const AutoDiff< SCAL> & y)
    {
        assert(x.Size() == y.Size());
        return x * Inv (y);
    }

    template< typename SCAL, typename SCAL2>
    inline AutoDiff< SCAL> operator/ (const AutoDiff< SCAL> & x, SCAL2 y)
    {
        return (1.0/y) * x;
    }

    template< typename SCAL, typename SCAL2>
    inline AutoDiff< SCAL> operator/ (SCAL2 x, const AutoDiff< SCAL> & y)
    {
        return x * Inv(y);
    }

    using std::sqrt;
    template< typename SCAL>
    inline AutoDiff< SCAL> sqrt (AutoDiff< SCAL> & x) throw()
    {
        AutoDiff< SCAL> res(x.Size(), 0);
        res.Value() = sqrt(x.Value());
        for (size_t j = 0; j < x.Size(); j++)  {
            res.DValue(j) = 0.5 / res.Value() * x.DValue(j);
        }

        return res;
    }

    // df(u)/dx  = exp(x) * du/dx
    // d^2 f(u) / dx^2 = exp(x) * (du/dx)^2 + exp(x) * d^2u /dx^2
    template < typename SCAL>
    inline AutoDiff< SCAL> exp (AutoDiff< SCAL> x)
    {
        AutoDiff< SCAL> res(x.Size(), 0);
        res.Value() = exp(x.Value());
        for (size_t k = 0; k < x.Size(); k++) {
            res.DValue(k) = x.DValue(k) * res.Value();
        }
        
        return res;
    }

    template < typename SCAL>
    inline AutoDiff< SCAL> log (AutoDiff< SCAL> x)
    {
        AutoDiff< SCAL> res(x.Size(), 0);
        res.Value() = log(x.Value());
        SCAL xinv = 1.0/x.Value();
        for (size_t k = 0; k < x.Size(); k++) {
            res.DValue(k) = x.DValue(k) * xinv;
        }
        
        return res;
    }

    using std::pow;
    template < typename SCAL>
    inline AutoDiff<SCAL> pow (AutoDiff<SCAL> x, AutoDiff<SCAL> y )
    {   
        assert(x.Size() == y.Size());
        return exp(log(x)*y);
    }

    using std::pow;
    template < typename SCAL>
    inline AutoDiff< SCAL> pow (AutoDiff< SCAL> x, SCAL y)
    {
        AutoDiff<SCAL> y_diff(x.Size(), y);
        
        return pow(x, y_diff);
    }

    template < typename SCAL>
    inline AutoDiff< SCAL> sin (AutoDiff< SCAL> x)
    {
        AutoDiff< SCAL> res(x.Size(), 0);
        SCAL s = sin(x.Value());
        SCAL c = cos(x.Value());
        
        res.Value() = s;
        for (size_t k = 0; k < x.Size(); k++) {
            res.DValue(k) = x.DValue(k) * c;
        }
        
        return res;
    }


    template < typename SCAL>
    inline AutoDiff< SCAL> cos (AutoDiff< SCAL> x)
    {
        AutoDiff< SCAL> res(x.Size(), 0);
        SCAL s = sin(x.Value());
        SCAL c = cos(x.Value());
        
        res.Value() = c;
        for (size_t k = 0; k < x.Size(); k++) {
            res.DValue(k) = -s * x.DValue(k);
        }
        
        return res;
    }

    template < typename SCAL>
    inline AutoDiff< SCAL> tan (AutoDiff< SCAL> x)
    { return sin(x) / cos(x); }
}

#endif