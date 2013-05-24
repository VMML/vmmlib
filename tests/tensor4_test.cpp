
#include "tensor4_test.hpp"

#include <vmmlib/tensor4.hpp>
#include <vmmlib/tucker4_tensor.hpp>
#include <vmmlib/t4_converter.hpp>
#include <vmmlib/tensor_mmapper.hpp>
#include <sstream>

#define TEST( x ) \
{ \
    ok = x; \
    global_ok &= ok; \
}

namespace vmml
{
    bool tensor4_test::run()
    {
        bool global_ok = true;
        bool ok = false;
        // indicates if failing the test produces only a warning
        bool fail_test = false;

        // define Tensor
        const size_t a = 3;
        const size_t b = 2;
        const size_t c = 4;
        const size_t d = 2;
        typedef long T;
        tensor4< a, b, c, d, T >  t4;

        // Test size
        TEST( t4.size() == a*b*c*d );
        if (ok)
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
            TEST( arrayptr[ index ] == 3);
            if (!ok) break;
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
            TEST( arrayptr[ index ] == 6);
            if (!ok) break;
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

        double check_value = 0;
        for ( size_t i4 = 0; i4 < d; ++i4 )
        {
            for ( size_t i3 = 0; i3 < c; ++i3 )
            {
                for( size_t i1 = 0; i1 < a; ++i1 )
                {
                    for( size_t i2 = 0; i2 < b && ok; ++i2 )
                    {
                        TEST(t4.at(i1, i2, i3, i4) == check_value);
                        if (!ok)
                        {
                            std::stringstream error;
                            error << "T4 with all values from 0 to 47:: " << std::endl << t4 << std::endl;
                            log_error( error.str() );
                        }
                        ++check_value;
                    }
                }
            }
        }

        if ( ok )
        {
            log( "tensor4 method fill_increasing()", true  );
        }

        ok = false;

        // Test indexing ( at and () )

        // test const versions as well
        const T c1 = t4(1,0,0,0);
        const T clast = t4.at(a-1,b-1,c-1,d-1);
        TEST( t4(0,0,0,0) == 0 &&
              t4.at(2, 1, 2, 1) == 1*a*b*c + 2*a*b + 1*a + 2 &&
              c1 == 2 && clast == T( t4.size( )) - 1 );
        if (ok)
        {
            log( "tensor4 indexing", true  );
        }
        else
        {
            std::stringstream error;
            error
            << "tensor4 indexing. t4 values should be: "
            << std::endl
            << 0 << "\t" << 1*a*b*c+2*a*b+1*a+2 << "\t" << 2 << "\t" << t4.size() - 1
            << std::endl
            << "t4 values are: " << std::endl
            << t4(0,0,0,0) << "\t" << t4.at(2, 1, 2, 1) << "\t" << c1 << "\t" << clast << std::endl;

            log_error( error.str(),  fail_test);
        }

       // Test equality
        tensor4<a, b, c, d, T> t4cp;
        t4.fill_increasing_values();
        t4cp.fill_increasing_values();
        tensor4<a, b, c, d, T> t4mod;
        t4mod.fill_increasing_values();
        t4mod(0,0,1,0) = 1;
        TEST( t4 == t4cp && ! (t4 != t4cp) && t4mod != t4 && !(t4mod == t4));
        if (ok)
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

        TEST( t4tempcp == t4 );
        if (ok)
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
            TEST((long)t4ptr[ index ] == t4longptr[ index ]);
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
        tensor4< a+2, b+2, c+2, d-1, T > t4resized(t4);
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

                        TEST(t4resized(i1, i2, i3, i4) == target);
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

        TEST(( (&t4assigned) != (&t4) ) && t4assigned == t4 );
        if (ok)
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

        TEST( !fail && max - min > (std::numeric_limits< T >::max)()/4); // assert at least a certain bandwith of values
        if (ok)
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

// FIXME rand() depends on stdlib version
//      ok = true;
//      tensor4< a, b, c, d, int >  t4_2;
//      t4_2.fill_random_signed(0);
//      int t4_test_data[] = {
//          -552808893, -183281951, -1044816132, -928209062, -250957408, 1058982018,
//          -33698213, -1007308382, 569808514, 929088271, -1005379225, 832964957,
//          196129103, -673590035, -45572427, 506060336, 715872987, -826667018,
//          426616161, 876625390, -310087809, -442380740, 296997003, -504719669,
//          -263079636, 1072593514, 89378518, -1046035170, -1053809277, 735507396,
//          756924037, 1023807655, -70843372, -651214229, -958621169, 766593553,
//          -764045132, 1043227471, 645667133, -679881061, 492627637, -514943,
//           -64720816, -868030333, 1011446114, 1018082584, -105720057, -235718411
//           };
//
//        size_t ind = 0;
//
//        for ( size_t i4 = 0; i4 < d; ++i4 )
//        {
//          for ( size_t i3 = 0; i3 < c; ++i3 )
//          {
//              for( size_t i1 = 0; i1 < a; ++i1 )
//              {
//                  for( size_t i2 = 0; i2 < b && ok; ++i2 )
//                  {
//                      TEST(t4_2.at(i1, i2, i3, i4) == t4_test_data[ind]);
//                      if (!ok)
//                      {
//                          std::stringstream error;
//                          error << "tensor4 random fill signed:: " << std::endl << t4_2 << std::endl;
//                          log_error( error.str() );
//                      }
//                      ++ind;
//                  }
//              }
//          }
//      }
//
//
//
//
//        if(min == 0 || max == 0 ) // at least one positive and one negative number
//        {
//            fail = true;
//        }
//
//        TEST( !fail && max - min > std::numeric_limits< T >::max()/2); // assert at least a certain bandwith of values
//        if (ok)
//      {
//          log( "tensor4 random fill signed", true  );
//      } else
//      {
//          std::stringstream error;
//          error
//          << "tensor4 random fill signed. t4 is: "
//          << std::endl
//          <<
//          t4_2
//          << std::endl;
//
//            if(fail) // if all values are the same or no negative/positive numbers, fail the test
//            {
//                log_error( error.str(),  fail_test);
//            }else
//            {
//                log_error( error.str(),  !fail_test);
//            }
//      }


      // Test casts
        ok = false;
        vmml::tensor4<a, b, c, d, ushort> t4same;
        t4same.fill_random();
        vmml::tensor4<a, b, c, d, T> t4control(t4same); // tested before so we can use the constructor
        t4.cast_from(t4same);
        TEST(t4 == t4control);
        if (ok)
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
        t4.fill_increasing_values();
        vmml::tensor4<a, b, c, d, ushort> t4diff;
        vmml::tensor4<a, b, c, d, T>    t4cast;
        t4diff.fill_increasing_values();
        t4cast.cast_from(t4diff);
        T val = 0;
        ok = true;


        for ( size_t i4 = 0; i4 < d; ++i4 )
        {
            for ( size_t i3 = 0; i3 < c; ++i3 )
            {
                for( size_t i1 = 0; i1 < a; ++i1 )
                {
                    for( size_t i2 = 0; i2 < b && ok; ++i2 )
                    {
                        //std::cout << "val is = " << t4cast.at(i1, i2, i3, i4) << "val should be = " << val << std::endl;
                        TEST(t4cast.at(i1, i2, i3, i4) == val);
                        if (!ok)
                        {
                            std::stringstream error;
                            error << "tensor4 cast_from different size test" << std::endl;
                            //error << "tensor4 cast_from different size: " << std::endl << t4 << std::endl;
                            log_error( error.str() );
                        }
                        ++val;
                    }
                }
            }
        }



        if(ok)
        {
            log( "tensor4 cast_from different size", true  );
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

        TEST( control == t4unsigned);
        if(ok)
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
            correctptr[ index ] = T( index+1 );
            TEST(calcptr[ index ] == T( index+1 ));
        }

        t4 += 1;

        if (ok) TEST(t4 == t4calc);
        if(ok)
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
            correctptr[ index ] = T( index-1 );
            TEST(calcptr[ index ] == T(index-1));
        }

        t4 -= 1;

        if(ok) TEST(t4 == t4calc);
        if (ok) {
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
            TEST(calcptr[ index ] == T( index*2 ));
        }

        t4 *= 2;

        if(ok) TEST(t4 == t4calc);
        if (ok)
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
            TEST(calcptr[ index ] == index/2);
        }

        t4 /= 2;

        if(ok) TEST(t4 == t4calc);
        if (ok)
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
            TEST(calcptr[ index ] == T( index*2 ));
        }

        t4 += t4offset;

        if(ok) TEST(t4 == t4calc);
        if (ok)
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
            TEST(calcptr[ index ] == 0);
        }

        t4 -= t4offset;

        if(ok) TEST(t4 == t4calc);
        if (ok)
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

        // FIXME
//        T* crrptr = t4origcrr.get_array_ptr();
//        T* addedptr = t4orig.get_array_ptr();
        //T* smallerptr = t4smaller.get_array_ptr();
//        ok = true;
//        T value = 0;
//        for(size_t index = 0; index < t4.size(); ++index)
//        {
//            crrptr[index] = ((index+1) % 3 == 0? index : index+value);
//            TEST( addedptr[index] ==  crrptr[index]);
//
//            if ((index +1) % 3 != 0)
//            {
//                ++value;
//            }
//        }
//
//        t4smaller.zero();
//        t4orig += t4smaller;
//
//        if(ok) TEST(t4orig == t4origcrr);
//        if (ok)
//        {
//            log( "tensor4 tensor different size addition +=", true  );
//        }else
//        {
//            std::stringstream error;
//          error
//          << "tensor4 tensor different size addition +=. t4 += is: "
//          << std::endl
//          << t4orig
//          << std::endl
//            << "t4 should be: " << std::endl
//            << t4origcrr
//            << std::endl;
//
//            log_error( error.str(),  fail_test);
//        }


       // Negation operators
        t4.fill_increasing_values();

        t4calc = -(t4);

        ok = true;
        for( size_t index = 0; index < t4calc.size(); ++index )
        {
            correctptr[index] = -ssize_t( index );
            TEST( ssize_t( calcptr[index] ) == -ssize_t( index ));
        }

        t4 = t4.negate();

        if(ok) TEST(t4 == t4calc);
        if (ok)
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

        //get_tensor 3 testing
        t4.fill_increasing_values();

        tensor3<a, b, c, T> t3 = t4.get_tensor3(0);
        tensor3<a, b, c, T> t3last = t4.get_tensor3(d-1);
        tensor3<a, b, c, T> t3correct;

        T offset = a*b*c*(d-1);
        ok = true;
        T* t3ptr = t3.get_array_ptr();
        T* t3lastptr = t3last.get_array_ptr();
        T* t3correctptr = t3correct.get_array_ptr();

        for(size_t index = 0; index < t3.size(); ++index)
        {
            t3correctptr[ index ] = index;
            TEST( t3ptr[index] == T( index ) &&
                  t3lastptr[index] == T( index + offset ));
        }

        if(ok)
        {
            log( "tensor4 get_tensor3()", true  );
        }
        else
        {
            std::stringstream error;
            error
            << "tensor4 get_tensor3(0). t3 is: "
            << std::endl
            << t3
            << std::endl
            << "t3 should be:"
            << std::endl
            <<  t3correct << std::endl;

            t3correct += offset;
            error
            << "tensor4 get_tensor3(n). t3 is: "
            << std::endl
            << t3last
            << std::endl
            << "t3 should be:"
            << std::endl
            <<  t3correct << std::endl;


            log_error( error.str(),  fail_test);
        }

      {
            //load mmap for tensor4

            //create test data
            std::string dir = ".";
            std::string filename = "mmap_testdata.raw";
            vmml::t4_converter<a,b,c,d,T> conv;

            t4.fill_increasing_values();
            conv.write_to_raw(t4, dir, filename);
            t4.fill_random(); // reset to make sure data is gone

            tensor_mmapper< tensor4<a,b,c,d,T>, vmml::t4_converter<a,b,c,d,T> > t4_mmap( dir, filename, false, conv );

            t4_mmap.get_tensor( t4 );

            tensor4< a,b,c,d,T > t4_check;
            t4_check.fill_increasing_values();

            TEST(t4_check == t4);
            if(ok)
            {
                log( "tensor4 load from memory mapped file", true  );
            }else
            {
                std::stringstream error;
                error
                << "tensor4 load from memory mapped file. t4 is: "
                << std::endl
                << t4
                << std::endl
                << "t4 should be:"
                << std::endl
                <<  t4_check << std::endl;

                log_error( error.str(),  fail_test);
            }


            remove("mmap_testdata.raw");
        }


       t4.fill_increasing_values();
        t3.zero();
        t3correct.zero();

        t4.average_I4(t3);
        ok = true;
        for (size_t index = 0; index < t3.size(); ++index)
        {
            t3correctptr[index] = ((d-1)*(a*b*c))/2+index;

            TEST(t3ptr[index] == t3correctptr[index]);
        }

        if (ok)
        {
            log( "tensor4 method average_I4()", true  );
        }else
        {
            std::stringstream error;
            error
            << "tensor4 method average_I4(). t3 is: "
            << std::endl
            << t3
            << std::endl
            << "t3 should be:"
            << std::endl
            << t3correct << std::endl;

            log_error( error.str(),  fail_test);

        }
        
        {
            tensor4< 2, 2, 2, 2, double > t4_data;
            double data[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
            t4_data.set( data, data + 16 );
            tucker4_tensor< 1, 1, 1, 1, 2, 2, 2, 2, double > tuck4_dec;
            typedef t4_hooi< 1, 1, 1, 1, 2, 2, 2, 2, double > hooi_type;
            tuck4_dec.tucker_als(t4_data, hooi_type::init_hosvd());
            tensor4< 2, 2, 2, 2, double > reco, reco_check;
            tuck4_dec.reconstruct( reco );
            
            double data_reco_check[] = { 2.5882, 2.8606, 3.1682, 3.5017, 3.9482, 4.3637, 4.8330, 5.3417, 7.7110, 8.5226, 9.4390, 10.4326, 11.7627, 13.0008, 14.3988, 15.9144 };
            reco_check.set( data_reco_check, data_reco_check + 16 );
            double precision = 1.0e-2;
            TEST( reco.equals(reco_check, precision) );
            if (ok)
            {
                log( "tensor4 Tucker reconstruction", true  );
            }else
            {
                std::stringstream error;
                error
                << "tensor4 Tucker reconstruction"
                << std::endl
                << "reco is:"
                << std::endl
                << reco
                << std::endl
                << "reco should be:"
                << std::endl
                << reco_check << std::endl;

                log_error( error.str(),  fail_test);
            }
        }

        return global_ok;
    }




} // namespace vmml
