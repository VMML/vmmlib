#include "tensor4_test.hpp"

#include <vmmlib/tensor4.hpp>
#include <sstream>

namespace vmml
{
	bool tensor4_test::run()
	{
		bool ok = false;
        // indicates if failing the test produces only a warning
        bool fail_test = false;

        // define Tensor
		const size_t a = 3;
        const size_t b = 2;
        const size_t c = 4;
        const size_t d = 4;
        typedef int T;
		tensor4< a, b, c, d, T >  t4;

        // Test size
        if ( t4.size() == a*b*c*d )
		{
			log( "tensor4 method size()", true  );
		} else
		{
			std::stringstream error;
			error
			<< "tensor4 method size(). size should be: "
			<< a*b*c*d
			<< std::endl
            << "size is: " << t4.size() << std::endl;
			log_error( error.str(),  fail_test);
		}

        // Test filling methods
        t4.fill(3);
        ok = true;
        T* arrayptr = t4.get_array_ptr();
        for(size_t index = 0; index < t4.size(); ++index)
        {
            if( arrayptr[ index ] != 3)
            {
                ok = false;
                break;
            }
        }
        if ( ok )
		{
			log( "tensor4 method fill( value )", true  );

		} else
		{
			std::stringstream error;
			error
			<< "tensor4 method fill( 3 ). t4 is: "
			<< std::endl
			<< t4
			<< std::endl;
			log_error( error.str(),  fail_test);
		}


        ok = true;
        t4 = 6;
        for(size_t index = 0; index < t4.size(); ++index)
        {
            if( arrayptr[ index ] != 6)
            {
                ok = false;
                break;
            }
        }
        if ( ok )
		{
			log( "tensor4 operator=( value )", true  );
		} else
		{
			std::stringstream error;
			error
			<< "tensor4 operator=( 6 ). t4 is: "
			<< std::endl
			<< t4
			<< std::endl;
			log_error( error.str(),  fail_test);
		}


        t4.fill_increasing_values();
        ok = true;
        for(size_t index = 0; index < t4.size(); ++index)
        {
            if( arrayptr[ index ] != index)
            {
                ok = false;
                break;
            }
        }

        if ( ok )
		{
			log( "tensor4 method fill_increasing()", true  );
		} else
		{
			std::stringstream error;
			error
			<< "tensor4 method fill_increasing(). t4 is: "
			<< std::endl
			<< t4
			<< std::endl;
			log_error( error.str(),  fail_test);
		}
        ok = false;

        // Test indexing ( at and () )

        // test const versions as well
        const int c1 = t4(1,0,0,0);
        const int clast = t4.at(a-1,b-1,c-1,d-1);
        if ( t4(0,0,0,0) == 0 && t4.at(2, 1, 2, 2) == 2*a*b*c + 2*a*b + 1*a + 2 && c1 == 1 && clast == t4.size() - 1 )
		{
			log( "tensor4 indexing", true  );
		} else
		{
			std::stringstream error;
			error
			<< "tensor4 indexing. t4 values should be: "
			<< std::endl
			<< 0 << "\t" << 2*a*b*c+2*a*b+a+2 << "\t" << 1 << "\t" << t4.size() - 1
			<< std::endl
            << "t4 values are: " << std::endl
            << t4(0,0,0,0) << "\t" << t4.at(2, 1, 2, 2) << "\t" << c1 << "\t" << clast << std::endl;

			log_error( error.str(),  fail_test);
		}

        // Test equality
        tensor4<a, b, c, d, T> t4cp;
        t4.fill_increasing_values();
        t4cp.fill_increasing_values();
        tensor4<a, b, c, d, T> t4mod;
        t4mod.fill_increasing_values();
        t4mod(0,0,1,0) = 1;
        if ( t4 == t4cp && ! (t4 != t4cp) && t4mod != t4 && !(t4mod == t4))
		{
			log( "tensor4 equality operators", true  );
		} else
		{
			std::stringstream error;
			error
			<< "tensor4 equality operators. Results should be: "
			<< std::endl
			<< true << true << true << true
			<< std::endl
            << "Results are:" << std::endl
            << (t4 == t4cp) << (! (t4 != t4cp)) << (t4mod != t4) << (! (t4mod == t4)) << std::endl;
			log_error( error.str(),  fail_test);
		}


        // Test copy ctors
        tensor4< a, b, c, d, T > t4tempcp(t4);

        if ( t4tempcp == t4 )
		{
			log( "tensor4 copy constructor", true  );
		} else
		{
			std::stringstream error;
			error
			<< "tensor4 copy constructor. t4 is: "
			<< std::endl
			<< t4
			<< std::endl
            << "tensor4 copy constructor. t4 copy is: "
            << std::endl
			<< t4tempcp
			<< std::endl;
			log_error( error.str(),  fail_test);
		}


        // cast constructor
        tensor4< a, b, c, d, long > t4long(t4);


        T* t4ptr = t4.get_array_ptr();
        long* t4longptr = t4long.get_array_ptr();
        ok = true;
        for(size_t index = 0; index < t4long.size() && ok; index++)
        {
            if((long)t4ptr[ index ] != t4longptr[ index ])
            {
                ok = false;
            }
        }

        if ( ok )
		{
			log( "tensor4 type converter constructor", true  );
		} else
		{
			std::stringstream error;
			error
			<< "tensor4 type converter constructor. t4 is: "
			<< std::endl
			<< t4
			<< std::endl
            << "tensor4 type converter constructor. t4 in long is: "
			<< std::endl
			<< t4long
			<< std::endl;
			log_error( error.str(),  fail_test);
		}

        ok = true;
        // resize constructor
        tensor4< a+2, b+2, c+2, d-2, T > t4resized(t4);
        T target = 0;
        for(size_t i4 = 0; i4 < t4resized.T3S; ++i4)
        {
            for(size_t i3 = 0; i3 < t4resized.SLICES; ++i3)
            {
                for(size_t i2 = 0; i2 < t4resized.COLS; ++i2)
                {
                    for(size_t i1 = 0; i1 < t4resized.ROWS; ++i1)
                    {
                        if(i1 < t4.ROWS && i2 < t4.COLS && i3 < t4.SLICES && i4 < t4.T3S )
                        {
                            target = t4(i1, i2, i3, i4);
                        }else
                        {
                            target = 0;
                        }

                        if(t4resized(i1, i2, i3, i4) != target)
                        {
                            ok = false;
                        }
                    }
                }
            }
        }

		if ( ok )
		{
			log( "tensor4 resize constructor", true  );
		} else
		{
			std::stringstream error;
			error
			<< "tensor4 resize constructor. t4 is: "
			<< std::endl
			<< t4
			<< std::endl
            << "tensor4 resize constructor. t4 resized is: "
			<< std::endl
			<< t4resized
			<< std::endl;
			log_error( error.str(),  fail_test);
		}
        ok = false;

        tensor4< a, b, c, d, T > t4assigned;
        t4assigned = t4;
        t4assigned = t4assigned; // also check self assignment

        if ( ( (&t4assigned) != (&t4) ) && t4assigned == t4 )
		{
			log( "tensor4 assignment operator=", true  );
		} else
		{
			std::stringstream error;
			error
			<< "tensor4 assignment operator=. t4 is: "
			<< std::endl
			<< t4 << " at " << &t4
			<< std::endl
            << "tensor4 assignment operator=. assigned t4 is: "
			<< std::endl
			<< t4assigned << " at " << &t4assigned
			<< std::endl;
			log_error( error.str(),  fail_test);
		}


		// fill random plausability tests
        t4.fill_random();
        t4ptr = t4.get_array_ptr();
        bool fail = true;
        T min = t4ptr[ 0 ];
        T max = t4ptr[ 0 ];
        for(size_t index = 0; index < t4.size(); index++)
        {
            if(t4ptr[ index ] < min)
            {
                min = t4ptr[ index ];
            }else if(t4ptr [ index ] > max)
            {
                max = t4ptr[ index ];
            }

            if(index > 0 && t4ptr[ index ] != t4ptr[ index-1 ]) // not all values can be identical
            {
                fail = false;
            }
        }

        if ( !fail && max - min > std::numeric_limits< T >::max()/4) // assert at least a certain bandwith of values
		{
			log( "tensor4 random fill", true  );
		} else
		{
			std::stringstream error;
			error
			<< "tensor4 random fill. t4 is: "
			<< std::endl
			<< t4
			<< std::endl;

            if(fail) // if all values are the same, fail the test
            {
                log_error( error.str(),  fail_test);
            }else
            {
                log_error( error.str(),  !fail_test);
            }
		}

        t4.fill_random_signed();
        t4ptr = t4.get_array_ptr();
        fail = true;
        min = 0;
        max = 0;
        for(size_t index = 0; index < t4.size(); index++)
        {
            if(t4ptr[ index ] < min)
            {
                min = t4ptr[ index ];
            }else if(t4ptr [ index ] > max)
            {
                max = t4ptr[ index ];
            }

            if(index > 0 && t4ptr[ index ] != t4ptr[ index-1 ]) // not all values can be identical
            {
                fail = false;
            }
        }

        if(min == 0 || max == 0 ) // at least one positive and one negative number
        {
            fail = true;
        }

        if ( !fail && max - min > std::numeric_limits< T >::max()/2) // assert at least a certain bandwith of values
		{
			log( "tensor4 random fill signed", true  );
		} else
		{
			std::stringstream error;
			error
			<< "tensor4 random fill signed. t4 is: "
			<< std::endl
			<< t4
			<< std::endl;

            if(fail) // if all values are the same or no negative/positive numbers, fail the test
            {
                log_error( error.str(),  fail_test);
            }else
            {
                log_error( error.str(),  !fail_test);
            }
		}

        // Test casts
        ok = false;
        vmml::tensor4<a, b, c, d, ushort> t4same;
        t4same.fill_random();
        vmml::tensor4<a, b, c, d, T> t4control(t4same); // tested before so we can use the constructor
        t4.cast_from(t4same);
        if(t4 == t4control)
        {
            log( "tensor4 cast_from same size", true  );
        }else
        {
            std::stringstream error;
			error
			<< "tensor4 cast_from same size. t4 is: "
			<< std::endl
			<< t4
			<< std::endl
            << "t4 should be: " << std::endl
            << t4control << std::endl;

            log_error( error.str(),  fail_test);
        }
        vmml::tensor4<4, 2, 4, 1, ushort> t4diff;
        vmml::tensor4<3, 2, 4, 4, T>    t4cast;
        t4diff.fill_increasing_values();
        t4cast.cast_from(t4diff);
        T val = 0;
        T* dat = t4cast.get_array_ptr();
        ok = true;
        for(size_t index = 0; index < 3*2*4 && ok; ++index)
        {
            if( (val+1) % 4 == 0)
            {
                ++val;
            }

            if(dat[index] != val)
            {
                ok = false;
            }

            ++val;
        }

        if(ok)
        {
            log( "tensor4 cast_from different size", true  );
        }else
        {
            std::stringstream error;
			error
			<< "tensor4 cast_from different size. t4 is: "
			<< std::endl
			<< t4
			<< std::endl
            << "t4 should be like: " << std::endl
            << t4diff << std::endl;

            log_error( error.str(),  fail_test);
        }


        vmml::tensor4<a, b, c, d, float> t4float;
        vmml::tensor4<a, b, c, d, unsigned short> t4unsigned;

        t4float.fill_increasing_values();
        float* datafloatptr = t4float.get_array_ptr();
        for(size_t index = 0; index < t4float.size(); ++index)
        {
            datafloatptr[index] = datafloatptr[index] + 0.45f;
        }

        t4unsigned.float_t_to_uint_t(t4float);
        vmml::tensor4<a, b, c, d, unsigned short> control;
        control.fill_increasing_values();

        if( control == t4unsigned)
        {
            log( "tensor4 float_t_to_uint_t", true  );
        }else
        {
            std::stringstream error;
			error
			<< "tensor4 float_t_to_uint_t. t4 is: "
			<< std::endl
			<< t4unsigned
			<< std::endl
            << "t4 should be: " << std::endl
            << control << std::endl;

            log_error( error.str(),  fail_test);
        }


        // Test scalar operators
        vmml::tensor4<a, b, c, d, T> t4calc;
        vmml::tensor4<a, b, c, d, T> t4correct;
        t4.fill_increasing_values();
        t4calc = t4 + 1;
        ok = true;
        T* calcptr = t4calc.get_array_ptr();
        T* correctptr = t4correct.get_array_ptr();
        for(size_t index = 0; index < t4calc.size(); ++index)
        {
            correctptr[ index ] = index+1;
            if(calcptr[ index ] != index+1)
            {
                ok = false;
                break;
            }
        }

        t4 += 1;

        if(ok && (t4 == t4calc))
        {
            log( "tensor4 scalar addition operators +/+=", true  );
        }else
        {
            std::stringstream error;
			error
			<< "tensor4 scalar addition operators +/+=. t4 += is: "
			<< std::endl
			<< t4
			<< std::endl
            << "t4 copy + is: "
			<< std::endl
			<< t4calc
			<< std::endl
            << "both should be: " << std::endl
            << t4correct << std::endl;

            log_error( error.str(),  fail_test);
        }


        t4.fill_increasing_values();
        t4calc = t4 - 1;
        ok = true;

        for(size_t index = 0; index < t4calc.size(); ++index)
        {
            correctptr[ index ] = index-1;
            if(calcptr[ index ] != index-1)
            {
                ok = false;
                break;
            }
        }

        t4 -= 1;

        if(ok && (t4 == t4calc))
        {
            log( "tensor4 scalar substraction operators -/-=", true  );
        }else
        {
            std::stringstream error;
			error
			<< "tensor4 scalar substraction operators -/-=. t4 -= is: "
			<< std::endl
			<< t4
			<< std::endl
            << "t4 copy - is: "
			<< std::endl
			<< t4calc
			<< std::endl
            << "both should be: " << std::endl
            << t4correct << std::endl;

            log_error( error.str(),  fail_test);
        }


        t4.fill_increasing_values();
        t4calc = t4 * 2;
        ok = true;

        for(size_t index = 0; index < t4calc.size(); ++index)
        {
            correctptr[ index ] = index*2;
            if(calcptr[ index ] != index*2)
            {
                ok = false;
                break;
            }
        }

        t4 *= 2;

        if(ok && (t4 == t4calc))
        {
            log( "tensor4 scalar multiplication operators */*=", true  );
        }else
        {
            std::stringstream error;
			error
			<< "tensor4 scalar multiplication operators */*=. t4 *= is: "
			<< std::endl
			<< t4
			<< std::endl
            << "t4 copy * is: "
			<< std::endl
			<< t4calc
			<< std::endl
            << "both should be: " << std::endl
            << t4correct << std::endl;

            log_error( error.str(),  fail_test);
        }

        t4.fill_increasing_values();
        t4calc = t4 / 2;
        ok = true;

        for(size_t index = 0; index < t4calc.size(); ++index)
        {
            correctptr[ index ] = index/2;
            if(calcptr[ index ] != index/2)
            {
                ok = false;
                break;
            }
        }

        t4 /= 2;

        if(ok && (t4 == t4calc))
        {
            log( "tensor4 scalar division operators / and /=", true  );
        }else
        {
            std::stringstream error;
			error
			<< "tensor4 scalar division operators / and /=. t4 /= is: "
			<< std::endl
			<< t4
			<< std::endl
            << "t4 copy / is: "
			<< std::endl
			<< t4calc
			<< std::endl
            << "both should be: " << std::endl
            << t4correct << std::endl;

            log_error( error.str(),  fail_test);
        }


        // Test Tensor scalar operators
        vmml::tensor4<a, b, c, d, T> t4offset;
        t4offset.fill_increasing_values();
        t4.fill_increasing_values();

        t4calc = t4 + t4offset;
        ok = true;

        for(size_t index = 0; index < t4calc.size(); ++index)
        {
            correctptr[ index ] = index*2;
            if(calcptr[ index ] != index*2)
            {
                ok = false;
                break;
            }
        }

        t4 += t4offset;

        if(ok && (t4 == t4calc))
        {
            log( "tensor4 tensor addition operators + and +=", true  );
        }else
        {
            std::stringstream error;
			error
			<< "tensor4 tensor addition operators + and +=. t4 += is: "
			<< std::endl
			<< t4
			<< std::endl
            << "t4 copy + is: "
			<< std::endl
			<< t4calc
			<< std::endl
            << "both should be: " << std::endl
            << t4correct << std::endl;

            log_error( error.str(),  fail_test);
        }

        t4offset.fill_increasing_values();
        t4.fill_increasing_values();

        t4calc = t4 - t4offset;
        ok = true;

        for(size_t index = 0; index < t4calc.size(); ++index)
        {
            correctptr[ index ] = 0;
            if(calcptr[ index ] != 0)
            {
                ok = false;
                break;
            }
        }

        t4 -= t4offset;

        if(ok && (t4 == t4calc))
        {
            log( "tensor4 tensor substraction operators - and -=", true  );
        }else
        {
            std::stringstream error;
			error
			<< "tensor4 tensor substraction operators - and -=. t4 -= is: "
			<< std::endl
			<< t4
			<< std::endl
            << "t4 copy - is: "
			<< std::endl
			<< t4calc
			<< std::endl
            << "both should be: " << std::endl
            << t4correct << std::endl;

            log_error( error.str(),  fail_test);
        }

        vmml::tensor4<2, 2, 4, 4, T> t4smaller;
        vmml::tensor4<3, 2, 4, 4, T> t4orig;
        vmml::tensor4<3, 2, 4, 4, T> t4origcrr;
        t4smaller.fill_increasing_values();
        t4orig.fill_increasing_values();
        t4orig += t4smaller;

        T* crrptr = t4origcrr.get_array_ptr();
        T* addedptr = t4orig.get_array_ptr();
        //T* smallerptr = t4smaller.get_array_ptr();

        ok = true;
        T value = 0;
        for(size_t index = 0; index < t4.size(); ++index)
        {
            crrptr[index] = ((index+1) % 3 == 0? index : index+value);
            if( addedptr[index] !=  crrptr[index])
            {
                ok = false;
            }

            if ((index +1) % 3 != 0)
            {
                ++value;
            }
        }

        t4smaller.zero();
        t4orig += t4smaller;

        if(ok && t4orig == t4origcrr)
        {
            log( "tensor4 tensor different size addition +=", true  );
        }else
        {
            std::stringstream error;
			error
			<< "tensor4 tensor different size addition +=. t4 += is: "
			<< std::endl
			<< t4orig
			<< std::endl
            << "t4 should be: " << std::endl
            << t4origcrr << std::endl;

            log_error( error.str(),  fail_test);
        }


        // Negation operators
        t4.fill_increasing_values();

        t4calc = -(t4);

        ok = true;
        for(size_t index = 0; index < t4calc.size(); ++ index)
        {
            correctptr[index] = -1 * (int) index;
            if(calcptr[index] != -1*(int)index )
            {
                ok = false;
            }
        }

        t4 = t4.negate();

        if(ok && (t4 == t4calc))
        {
            log( "tensor4 tensor negation", true  );
        }else
        {
            std::stringstream error;
			error
			<< "tensor4 tensor negation. t4.negate() is: "
			<< std::endl
			<< t4
			<< std::endl
            << "-(t4) is: "
			<< std::endl
			<< t4calc
			<< std::endl
            << "both should be: " << std::endl
            << t4correct << std::endl;

            log_error( error.str(),  fail_test);
        }

		return ok;
	}



} // namespace vmml

