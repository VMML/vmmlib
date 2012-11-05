/* 
 * VMMLib - Tensor Classes
 *  
 * @author Susanne Suter
 * @author Rafa Ballester
 *
 * The higher-order orthogonal iteration (HOOI) is also known as Tucker-ALS (Tuck-ALS)
 * The t3_hooi implements a HOOI for a third-order tensor
 * references:
 * - Tucker, 1966: Some mathematical notes on three-mode factor analysis, Psychometrika.
 * - De Lathauwer, De Moor, Vandewalle, 2000a: A multilinear singular value decomposition, SIAM J. Matrix Anal. Appl.
 * - De Lathauwer, De Moor, Vandewalle, 2000b: On the Best rank-1 and Rank-(R_1, R_2, ..., R_N) Approximation and Applications of Higher-Order Tensors, SIAM J. Matrix Anal. Appl.
 * - Kolda & Bader, 2009: Tensor Decompositions and Applications, SIAM Review.
 * - Bader & Kolda, 2006: Algorithm 862: Matlab tensor classes for fast algorithm prototyping. ACM Transactions on Mathematical Software.
 * 
 */

#ifndef __VMML__T3_HOOI__HPP__
#define __VMML__T3_HOOI__HPP__


#include <vmmlib/t3_hosvd.hpp>
#include <vmmlib/t3_ttm.hpp>
#include <vmmlib/matrix_pseudoinverse.hpp>
#include <vmmlib/tensor_stats.hpp>

namespace vmml {

    template< size_t R1, size_t R2, size_t R3, size_t I1, size_t I2, size_t I3, typename T = float >
            class t3_hooi {
    public:

        typedef tensor3< I1, I2, I3, T > t3_type;

        typedef tensor3< R1, R2, R3, T > t3_core_type;

        typedef matrix< I1, R1, T > u1_type;
        typedef matrix< I2, R2, T > u2_type;
        typedef matrix< I3, R3, T > u3_type;

        typedef matrix< R1, I1, T > u1_t_type;
        typedef matrix< R2, I2, T > u2_t_type;
        typedef matrix< R3, I3, T > u3_t_type;

        /*	higher-order orthogonal iteration (HOOI) is a truncated HOSVD decompositions, i.e., the HOSVD components are of lower-ranks. An optimal rank-reduction is 
         performed with an alternating least-squares (ALS) algorithm, which minimizes the error between the approximated and orignal tensor based on the Frobenius norm
         see: De Lathauwer et al, 2000b; On the best rank-1 and rank-(RRR) approximation of higher-order tensors.
         the HOOI can be computed based on (a) n-mode PCA, i.e., an eigenvalue decomposition on the covariance matrix of every mode's matriciziation, and 
         (b) by performing a 2D SVD on the matricization of every mode. Matrix matricization means that a tensor I1xI2xI3 is unfolded/sliced into one matrix
         with the dimensions I1xI2I3, which corresponds to a matrizitation alonge mode I1.
         */
        template< typename T_init>
        static tensor_stats als(const t3_type& data_, u1_type& u1_, u2_type& u2_, u3_type& u3_, t3_core_type& core_, T_init init, const double& max_f_norm_ = 0.0, const size_t max_iterations = 3, const float tolerance = 1e-04);

        //core not needed
        template< typename T_init>
        static tensor_stats als(const t3_type& data_, u1_type& u1_, u2_type& u2_, u3_type& u3_, T_init init, const double& max_f_norm_ = 0.0, const size_t max_iterations = 3, const float tolerance = 1e-04);

        /* derive core
         implemented according to core = data x_1 U1_pinv x_2 U2_pinv x_3 U3_pinv, 
         where x_1 ... x_3 are n-mode products and U1_pinv ... U3_pinv are inverted basis matrices
         the inversion is done with a matrix pseudoinverse computation
         */
        static void derive_core(const t3_type& data_, const u1_type& u1_, const u2_type& u2_, const u3_type& u3_, t3_core_type& core_);
        //faster: but only if basis matrices are orthogonal
        static void derive_core_orthogonal_bases(const t3_type& data_, const u1_type& u1_, const u2_type& u2_, const u3_type& u3_, t3_core_type& core_);

        // init functors 

        struct init_hosvd {

            inline void operator()(const t3_type& data_, u1_type& u1_, u2_type& u2_, u3_type & u3_) {
                t3_hosvd< R1, R2, R3, I1, I2, I3, T >::apply_mode1(data_, u1_);
                t3_hosvd< R1, R2, R3, I1, I2, I3, T >::apply_mode2(data_, u2_);
                t3_hosvd< R1, R2, R3, I1, I2, I3, T >::apply_mode3(data_, u3_);
            }
        };

        struct init_random {

            inline void operator()(const t3_type& data_, u1_type& u1_, u2_type& u2_, u3_type & u3_) {
                srand(time(NULL));
                u1_.set_random();
                u2_.set_random();
                u3_.set_random();

                u1_ /= u1_.frobenius_norm();
                u2_ /= u2_.frobenius_norm();
                u3_ /= u3_.frobenius_norm();
            }
        };

        struct init_dct {

            inline void operator()(const t3_type& data_, u1_type& u1_, u2_type& u2_, u3_type & u3_) {
                u2_.set_dct();
                u3_.set_dct();
            }
        };

    protected:


        static void optimize_mode1(const t3_type& data_, const u2_type& u2_, const u3_type& u3_,
                tensor3< I1, R2, R3, T >& projection_,
                tensor3< I1, R2, I3, T >& tmp_);
        static void optimize_mode2(const t3_type& data_, const u1_type& u1_, const u3_type& u3_,
                tensor3< R1, I2, R3, T >& projection_,
                tensor3< R1, I2, I3, T >& tmp_);
        static void optimize_mode3(const t3_type& data_, const u1_type& u1_, const u2_type& u2_,
                tensor3< R1, R2, I3, T >& projection_,
                tensor3< R1, I2, I3, T >& tmp_);


    }; //end class t3_hooi

#define VMML_TEMPLATE_STRING        template< size_t R1, size_t R2, size_t R3, size_t I1, size_t I2, size_t I3, typename T >
#define VMML_TEMPLATE_CLASSNAME     t3_hooi< R1, R2, R3, I1, I2, I3, T >

    VMML_TEMPLATE_STRING
    template< typename T_init>
    tensor_stats
    VMML_TEMPLATE_CLASSNAME::als(const t3_type& data_,
            u1_type& u1_, u2_type& u2_, u3_type& u3_,
            T_init init, const double& max_f_norm_, const size_t max_iterations, const float tolerance) {
        t3_core_type core;
        core.zero();
        return als(data_, u1_, u2_, u3_, core, init, max_f_norm_, max_iterations, tolerance);
    }

    VMML_TEMPLATE_STRING
    template< typename T_init>
    tensor_stats
    VMML_TEMPLATE_CLASSNAME::als(const t3_type& data_,
            u1_type& u1_, u2_type& u2_, u3_type& u3_,
            t3_core_type& core_,
            T_init init,
            const double& max_f_norm_, const size_t max_iterations_, const float tolerance_) {
        tensor_stats result;

        //intialize basis matrices
        init(data_, u1_, u2_, u3_);

        core_.zero();
        double max_f_norm, f_norm, fit, fitchange, fitold, normresidual;
        if (tolerance_ > 0) {
            max_f_norm = max_f_norm_;

            if (max_f_norm <= 0.0) {
                max_f_norm = data_.frobenius_norm();
            }
            fit = 0;
            //removed to save computation
            /*if ( (max_f_norm != 0) && (max_f_norm > f_norm) ) 
            {
                    fit = 1 - (normresidual / max_f_norm);
            } else { 
                    fit = 1;
            }*/
            fitchange = 1;
            fitold = fit;
            normresidual = 0;
        }

        tensor3< I1, R2, R3, T > projection1;
        tensor3< R1, I2, R3, T > projection2;
        tensor3< R1, R2, I3, T > projection3;

        tensor3< I1, R2, I3, T > tmp1;
        tensor3< R1, I2, I3, T > tmp2;

#if TUCKER_LOG
        std::cout << "HOOI ALS (for tensor3) " << std::endl
                << "initial fit: " << fit << ", "
                << "frobenius norm original: " << max_f_norm << std::endl;
#endif	
        size_t i = 0;
        while (i < max_iterations_ && (tolerance_ == -1 || fitchange >= tolerance_)) { //do until converges
            fitold = fit;

            //optimize modes
            optimize_mode1(data_, u2_, u3_, projection1, tmp1);
            t3_hosvd< R1, R2, R3, I1, R2, R3, T >::apply_mode1(projection1, u1_);

            optimize_mode2(data_, u1_, u3_, projection2, tmp2);
            t3_hosvd< R1, R2, R3, R1, I2, R3, T >::apply_mode2(projection2, u2_);

            optimize_mode3(data_, u1_, u2_, projection3, tmp2);
            t3_hosvd< R1, R2, R3, R1, R2, I3, T >::apply_mode3(projection3, u3_);

            t3_ttm::multiply_horizontal_bwd(projection3, transpose(u3_), core_);

            if (tolerance_ > 0) {
                f_norm = core_.frobenius_norm();
                normresidual = sqrt(max_f_norm * max_f_norm - f_norm * f_norm);
                fit = 1 - (normresidual / max_f_norm);
                fitchange = fabs(fitold - fit);
#if TUCKER_LOG
                std::cout << "iteration '" << i << "', fit: " << fit
                        << ", fitdelta: " << fitchange
                        << ", frobenius norm of core: " << f_norm << std::endl;
#endif
            }
            ++i;
        }
        result.set_n_iterations(i);
        return result;
    }

    VMML_TEMPLATE_STRING
    void
    VMML_TEMPLATE_CLASSNAME::optimize_mode1(const t3_type& data_, const u2_type& u2_, const u3_type& u3_,
            tensor3< I1, R2, R3, T >& projection_,
            tensor3< I1, R2, I3, T >& tmp_) {
        u2_t_type* u2_inv = new u2_t_type;
        u3_t_type* u3_inv = new u3_t_type;
        u2_.transpose_to(*u2_inv);
        u3_.transpose_to(*u3_inv);

#if 1
        //backward cyclic matricization/unfolding (after Lathauwer et al., 2000a)
        t3_ttm::multiply_frontal_bwd(data_, *u2_inv, tmp_);
        t3_ttm::multiply_horizontal_bwd(tmp_, *u3_inv, projection_);
#else
        //forward cyclic matricization/unfolding (after Kiers, 2000) -> memory optimized
        t3_ttm::multiply_horizontal_fwd(data_, *u2_inv, tmp_);
        t3_ttm::multiply_lateral_fwd(tmp_, *u3_inv, projection_);
#endif

        delete u2_inv;
        delete u3_inv;
    }

    VMML_TEMPLATE_STRING
    void
    VMML_TEMPLATE_CLASSNAME::optimize_mode2(const t3_type& data_, const u1_type& u1_, const u3_type& u3_,
            tensor3< R1, I2, R3, T >& projection_,
            tensor3< R1, I2, I3, T >& tmp_) {
        u1_t_type* u1_inv = new u1_t_type();
        u3_t_type* u3_inv = new u3_t_type();
        u1_.transpose_to(*u1_inv);
        u3_.transpose_to(*u3_inv);


#if 0
        //backward cyclic matricization (after Lathauwer et al., 2000a)
        t3_ttm::multiply_lateral_bwd(data_, *u1_inv, tmp_);
        t3_ttm::multiply_horizontal_bwd(tmp_, *u3_inv, projection_);
#else
        //forward cyclic matricization/unfolding (after Kiers, 2000) -> memory optimized
        t3_ttm::multiply_frontal_fwd(data_, *u1_inv, tmp_);
        t3_ttm::multiply_lateral_fwd(tmp_, *u3_inv, projection_);
#endif

        delete u1_inv;
        delete u3_inv;
    }

    VMML_TEMPLATE_STRING
    void
    VMML_TEMPLATE_CLASSNAME::optimize_mode3(const t3_type& data_, const u1_type& u1_, const u2_type& u2_,
            tensor3< R1, R2, I3, T >& projection_,
            tensor3< R1, I2, I3, T >& tmp_) {
        u1_t_type* u1_inv = new u1_t_type;
        u2_t_type* u2_inv = new u2_t_type;
        u1_.transpose_to(*u1_inv);
        u2_.transpose_to(*u2_inv);

#if 0
        //backward cyclic matricization (after Lathauwer et al., 2000a)
        t3_ttm::multiply_lateral_bwd(data_, *u1_inv, tmp_);
        t3_ttm::multiply_frontal_bwd(tmp_, *u2_inv, projection_);
#else
        //forward cyclic matricization/unfolding (after Kiers, 2000) -> memory optimized
        t3_ttm::multiply_frontal_fwd(data_, *u1_inv, tmp_);
        t3_ttm::multiply_horizontal_fwd(tmp_, *u2_inv, projection_);
#endif

        delete u1_inv;
        delete u2_inv;
    }

    VMML_TEMPLATE_STRING
    void
    VMML_TEMPLATE_CLASSNAME::derive_core_orthogonal_bases(const t3_type& data_, const u1_type& u1_, const u2_type& u2_, const u3_type& u3_, t3_core_type& core_) {
        u1_t_type* u1_inv = new u1_t_type;
        u2_t_type* u2_inv = new u2_t_type;
        u3_t_type* u3_inv = new u3_t_type;

        u1_.transpose_to(*u1_inv);
        u2_.transpose_to(*u2_inv);
        u3_.transpose_to(*u3_inv);

        t3_ttm::full_tensor3_matrix_multiplication(data_, *u1_inv, *u2_inv, *u3_inv, core_);

        delete u1_inv;
        delete u2_inv;
        delete u3_inv;
    }

    VMML_TEMPLATE_STRING
    void
    VMML_TEMPLATE_CLASSNAME::derive_core(const t3_type& data_, const u1_type& u1_, const u2_type& u2_, const u3_type& u3_, t3_core_type& core_) {

#if 1
        //faster approach
        //compute pseudo inverse for matrices u1-u3
        u1_type* u1_pinv_t = new u1_type;
        u2_type* u2_pinv_t = new u2_type;
        u3_type* u3_pinv_t = new u3_type;

        compute_pseudoinverse< u1_type > compute_pinv_u1;
        compute_pinv_u1(u1_, *u1_pinv_t);
        compute_pseudoinverse< u2_type > compute_pinv_u2;
        compute_pinv_u2(u2_, *u2_pinv_t);
        compute_pseudoinverse< u3_type > compute_pinv_u3;
        compute_pinv_u3(u3_, *u3_pinv_t);

        u1_t_type* u1_pinv = new u1_t_type;
        u2_t_type* u2_pinv = new u2_t_type;
        u3_t_type* u3_pinv = new u3_t_type;

        u1_pinv_t->transpose_to(*u1_pinv);
        u2_pinv_t->transpose_to(*u2_pinv);
        u3_pinv_t->transpose_to(*u3_pinv);

        t3_ttm::full_tensor3_matrix_multiplication(data_, *u1_pinv, *u2_pinv, *u3_pinv, core_);

        delete u1_pinv;
        delete u2_pinv;
        delete u3_pinv;
        delete u1_pinv_t;
        delete u2_pinv_t;
        delete u3_pinv_t;

#else
        //previous version of compute core	
        for (size_t r3 = 0; r3 < R3; ++r3) {
            for (size_t r1 = 0; r1 < R1; ++r1) {
                for (size_t r2 = 0; r2 < R2; ++r2) {
                    float_t sum_i1_i2_i3 = 0.0;
                    for (size_t i3 = 0; i3 < I3; ++i3) {
                        for (size_t i1 = 0; i1 < I1; ++i1) {
                            for (size_t i2 = 0; i2 < I2; ++i2) {
                                sum_i1_i2_i3 += u1_.at(i1, r1) * u2_.at(i2, r2) * u3_.at(i3, r3) * T(data_.at(i1, i2, i3));
                            }
                        }
                    }
                    core_.at(r1, r2, r3) = sum_i1_i2_i3;
                }
            }
        }

#endif
    }


#undef VMML_TEMPLATE_STRING
#undef VMML_TEMPLATE_CLASSNAME


}//end vmml namespace

#endif

